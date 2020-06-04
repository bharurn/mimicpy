# This is part of MiMiCPy

"""

This module contains the handles to prepare the
MM topology, MPT and QM region/CPMD script

"""

from .selector import Selector
from .base import BaseHandle
from .._global import _Global as _global
from . import _qmhelper
from ..parsers.mpt import Reader as MPTReader, write as mptwrite
from ..parser._mpt_writer import atomtypes
from ..utils.constants import hartree_to_ps, bohr_rad
from ..scripts import cpmd, mdp
from ..parsers import parse_pdb

class MM(BaseHandle):
    """
    Prepare MM topology, by writing MPT file, also supports non-std residues
    Inherits from .core.base.BaseHandle
    
    """
    
    def __init__(self, topol_dir='prepareMM', status=None):
        """Class constructor"""
        
        super().__init__(status) # call BaseHandle constructor to init _status dict
        
        if self._status['prepMM'] != '':
            self.dir = self._status['prepMM']
        else:
            self.dir = topol_dir
            
        _global.logger.write('debug', f'Set handle directory as {self.dir}..')
        
        self.toYaml()
        
        self.nonstd_atm_types = {}
    
    def addLig(self, pdb, itp, *resnames, buff=1000):
        pdb_df = parse_pdb.parseFile(pdb, buff)
        
        elems = [pdb_df['element'][pdb_df['resName']==name] for name in resnames]
        
        f = _global.host.open(itp, 'rb')
        
        atm_types = atomtypes(f, buff, True)
        f.close()
        
        self.nonstd_atm_types.update( dict(zip(atm_types, elems)) )
        
    
    def getMPT(self, mpt=None, buff=1000, guess_elems=False):
        """Get the MPT topology, used in prepare.QM"""
        
        # generate proprocessed topology for now, TO DO: changed writer mpt to read from .top
        self.grompp(mdp.MDP.defaultGenion(), self.ions, gro = self.conf2, pp = self.preproc, dirc=self.dir) 
        
        if mpt == None: mpt = f"{self.dir}/{self.mpt}" # if no mpt file was passed, use default value
        else: mpt = f"{self.dir}/{mpt}"
        
        mptwrite(self.preproc, mpt, self.nonstd_atm_types, buff, guess_elems)
        
        self.toYaml()
        
class QM(BaseHandle):
    """
    Prepare QM input script, by reading MPT file
    Inherits from .core.base.BaseHandle
    
    """
    
    def __init__(self, status=None, selector=Selector(), mpt=None, gro=None):
        """Class constructor"""
        
        super().__init__(status) # call BaseHandle.__init__() to init status dict
        
        
        self.mpt = MPTReader(self.getcurrentNone(mpt, 'mpt'))
        self.gro = self.getcurrentNone(gro, 'gro')
        self.selector = selector
        
        # load mpt and gro into selector, returns box size in gro
        self._mm_box = self.selector.load(self.mpt, self.gro)
        
        # init scripts and paths/files
        self.inp = cpmd.Input()
        self.qmatoms = None
        if self._status['prepQM'] == '':
            self.dir = 'prepareQM'
        else:
            self.dir = self._status['prepQM']
        self.index = 'index.ndx'
        self.preprc = 'processed.top'
        self.mimic = 'mimic'
    
    def add(self, selection=None, link=False):
        """Add dataframe to QM region"""
        
        qdf = _qmhelper._cleanqdf( self.selector.select(selection) )
        
        # add a new column link which is integer of link argument
        qdf.insert(2, 'link', [int(link)]*len(qdf))
        
        # add qdf to self.qmatoms, append if already exists
        if self.qmatoms is None:
            self.qmatoms = qdf
        else:
            self.qmatoms = self.qmatoms.append(qdf)
    
    def delete(self, selection=None):
        """Delete from QM region using selection langauage"""
        qdf = QM._cleanqdf( self.selector.select(selection) )
        # drop selection, ignore errors to ignore extra residues selected that are not present in self.qmatoms
        self.qmatoms = self.qmatoms.drop(qdf.index, errors='ignore')
    
    def clear(self):
        """Empty the QM region to start over"""
        self.qmatoms = None
    
    def getInp(self, mdp=mdp.MDP.defaultMiMiC()):
        """
        Create the QM region from the atoms added
        Steps done:
            Set integrator to mimic, and QMMM_grps to QMatoms in mdp script
            Writes index files of atoms in QM region with name QMatoms
            Run gmx grompp to get mimic.tpr
            Sort self.qmatoms by link and elements
            Read element symbols and all charge info from self.qmatoms
            Add MIMIC, SYSTEM section of CPMD script
            Add all atoms to CPMD script
        """
        dirc = self.dir
        self.setcurrent(key='prepQM')
        
        _global.logger.write('debug2', "Changing Gromacs integrator to MiMiC..")
        mdp.integrator = 'mimic'
        
        _global.logger.write('debug', f"Writing atoms in QM region to {self.index}..")
        mdp.QMMM_grps = 'QMatoms'
        _global.host.write(_qmhelper.index(self.qmatoms.index, mdp.QMMM_grps), f'{dirc}/{self.index}')
        
        _global.logger.write('info', "Generating Gromacs tpr file for MiMiC run..")
        self.grompp(mdp, self.mimic, gro=self.gro, n=self.index, dirc=dirc)
        
        # sort by link column first, then element symbol
        # ensures that all link atoms are last, and all elements are bunched together
        # index is also reset, for getOverlaps_Atoms()
        sorted_qm = self.qmatoms.sort_values(by=['link', 'element']).reset_index()
        
        _global.logger.write('info', "Creating CPMD input script..")
        inp=cpmd.Input()
        
        inp.mimic = cpmd.Section()
        inp.mimic.paths = "1\n---" #path will be set by simulate.MiMiC.run()
        inp.mimic.box = '  '.join([str(s/bohr_rad) for s in self._mm_box])
        
        inp = _qmhelper.getOverlaps_Atoms(sorted_qm, inp) # get the overlaps and add atoms section of inp
        
        inp.mimic.long_range__coupling = ''
        inp.mimic.fragment__sorting__atom_wise__update = 100
        inp.mimic.cutoff__distance = 20.0
        inp.mimic.multipole__order = 3
        
        q = sum(self.qmatoms['charge'])
        if not round(q, 2).is_number(): 
            _global.logger.write('warning', (f'Total charge of QM region ({q}) not an integer up to 2 decimal places.'
                                            ' \nRounding to integer anyways..'))
        
        inp.system.charge = round(q) # system section already created in getQverlap_Atoms()
        
        # set cpmd section
        if not inp.checkSection('cpmd'): inp.cpmd = cpmd.Section()
        inp.cpmd.mimic = ''
        inp.cpmd.parallel__constraints = ''
        
        # set no of steps from mdp file
        if not mdp.hasparam('nsteps'): mdp.nsteps = 1000 # default value, if not present in mdp
        inp.cpmd.maxsteps = mdp.nsteps
        
        inp.dft = cpmd.Section()
        inp.dft.functional__blyp = '' # update to new XC_DRIVER code
        
        # set timestep from mdp file
        if not mdp.hasparam('dt'): mdp.dt = 0.0001 # default value, if not present in mdp
        inp.cpmd.timestep = round(mdp.dt/hartree_to_ps) # convert to atomic units and round off
        
        self.inp = inp
        
        self.toYaml()
        
        return inp
