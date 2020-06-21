# This is part of MiMiCPy

"""

This module contains the handles to prepare the
MM topology, MPT and QM region/CPMD script

"""

from .selector import Selector
from .base import BaseHandle
from .._global import _Global as _global
from . import _qmhelper
from ..parsers.mpt import MPT
from ..parsers.top_reader import ITPParser
from ..utils.constants import hartree_to_ps, bohr_rad
from ..scripts import cpmd, mdp
from ..parsers import pdb as parse_pdb
from ..utils.errors import MiMiCPyError

class MM(BaseHandle):
    """
    Prepare MM topology, by writing MPT file, also supports non-std residues
    Inherits from .core.base.BaseHandle
    
    """
    
    def __init__(self, topol_dir='', status=None, savestatus=True):
        """Class constructor"""
        
        super().__init__(status) # call BaseHandle constructor to init _status dict
        
        self.savestatus = savestatus
        
        if self._status['prepMM'] != '':
            self.dir = self._status['prepMM']
        else:
            self.dir = topol_dir
            
        if self.dir:
            _global.logger.write('debug', f'Set handle directory as {self.dir}..')
        else:
            _global.logger.write('debug', f'Set handle directory as current directory..')
        
        self.toYaml()
        
        self.nonstd_atm_types = {}
        self.conf1 = 'solv.gro'
        self.conf2 = 'ions.gro'
        self.ions = "ions"
        
        self._ion_kwargs = {'pname': 'NA', 'nname': 'CL', 'neutral': ''} # parameters to pass to gmx genion
        self._solavte_kwargs = {'cs': 'spc216.gro'} # parameters to pass to gmx solvate
    
    def addLig(self, pdb, itp, *resnames, buff=1000):
        pdb_df = parse_pdb.parseFile(pdb, buff)
        
        # get chain info
        chains = pdb_df['chainID'].unique().tolist()
        chains.remove(' ') # remove elements wtih no chain info from lisy
        
        # either select only first chain of residue or all residues with no chain info
        pdb_df = pdb_df[(pdb_df['chainID'] == chains[0]) | (pdb_df['chainID'] == ' ')]
        
        # get elems
        elems = [a for name in resnames for a in pdb_df['element'][pdb_df['resName']==name].to_list()]
        
        # get atom types
        ITPParser.clear()
        
        # fake class to satisfy atomtypes dict
        class defdict:
            def __contains__(self, item): return True
            def __getitem__(self,key): return ' '
        
        itp_parser = ITPParser(resnames, defdict(), buff, False)
        itp_parser.parse(itp)
        atm_types = [a for df in itp_parser.dfs for a in df['type'].to_list()]
        
        # assert len(atm_types) == len(elems) before zipping
        self.nonstd_atm_types.update( dict(zip(atm_types, elems)) )
        
    
    def getMPT(self, topol=None, mpt=None, buff=1000, guess_elems=False, toFile=True):
        """Get the MPT topology, used in prepare.QM"""
        
        top = self.getcurrentNone(topol, 'top')
        
        if not toFile: return MPT.fromTop(top, self.nonstd_atm_types, buff, guess_elems)
        else: mpt_handle = MPT.fromTop(top, self.nonstd_atm_types, buff, guess_elems, mode='w')
        
        if mpt == None: mpt = _global.host.join(self.dir, 'topol.mpt') # if no mpt file was passed, use default value
        else: mpt = _global.host.join(self.dir,  mpt)
        
        mpt_handle.write(mpt)
        
        self.toYaml()
    
    def genionParams(self, **kwargs): self._ion_kwargs = kwargs
    def solvateParams(self, **kwargs): self._solavte_kwargs = kwargs
    
    def makeBox(self, gro=None, genion_mdp=mdp.MDP.defaultGenion()):
        """Run gmx solvate, gmx grompp and gmx genion in that order"""
        
        _global.logger.write('info', "Solvating box..")
        self.gmx('solvate', cp = self.getcurrentNone(gro, 'gro'), o = self.conf1, p = self.topol, dirc=self.dir, **self._solavte_kwargs)
        
        _global.logger.write('info', "Adding ions to box..")
        self.grompp(genion_mdp, self.ions, gro = self.conf1, dirc=self.dir)
        
        # send SOL to stdin, so gromacs replaces some SOL molecules with ions 
        self.gmx('genion', s = self.ions_tpr, o = self.conf2, p = self.topol, dirc=self.dir, **self._ion_kwargs, stdin="SOL")
        
        _global.logger.write('info', 'Simulation box prepared..')
        
        self.saveToYaml()
        
class QM(BaseHandle):
    """
    Prepare QM input script, by reading MPT file
    Inherits from .core.base.BaseHandle
    
    """
    
    def __init__(self, status=None, selector=Selector(), mpt=None, gro=None):
        """Class constructor"""
        
        super().__init__(status) # call BaseHandle.__init__() to init status dict
        
        if isinstance(mpt, MPT): self.mpt = mpt
        else: self.mpt = MPT.fromFile(self.getcurrentNone(mpt, 'mpt'))
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
        self.QMMM_grps = 'QMatoms'
    
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
        # drop selection, ignore pandas errors to disregard extra residues selected that are not present in self.qmatoms
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
            Add MIMIC, SYSTEM, DFT sections of CPMD script
            Add all atoms to CPMD script
        """
        if self.qmatoms == None:
            raise MiMiCPyError("No QM atoms have been selected")
        
        self.mpt.close() # clear all_data from memory, in case gc doesn't work
        ####
        #Write ndx, tpr file for gromacs
        #output only ndx if mdp is None
        ####
        dirc = self.dir
        self.setcurrent(key='prepQM')
        
        # write index file
        _global.host.write(_qmhelper.index(self.qmatoms.index, self.QMMM_grps), f'{dirc}/{self.index}')
        
        # default vals
        nsteps = 1000
        dt = 0.0001
            
        if mdp != None:
            _global.logger.write('debug2', "Changing Gromacs integrator to MiMiC..")
            mdp.integrator = 'mimic'
        
            _global.logger.write('debug', f"Writing atoms in QM region to {self.index}..")
            mdp.QMMM_grps = self.QMMM_grps
            
            _global.logger.write('info', "Generating Gromacs .tpr file for MiMiC run..")
            
            # set no of steps in mdp file
            if not mdp.hasparam('nsteps'): mdp.nsteps = nsteps # default value, if not present in mdp
            else: nsteps = mdp.nsteps
            
            # set timestep in mdp file
            if not mdp.hasparam('dt'): mdp.dt = dt # default value, if not present in mdp
            else: dt = mdp.dt
            
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
        if not round(q, 2).is_integer(): 
            _global.logger.write('warning', (f'Total charge of QM region ({q}) not an integer up to 2 decimal places.'
                                            ' \nRounding to integer anyways..'))
        
        inp.system.charge = round(q) # system section already created in getQverlap_Atoms()
        
        # set cpmd section
        if not inp.checkSection('cpmd'): inp.cpmd = cpmd.Section()
        inp.cpmd.mimic = ''
        inp.cpmd.parallel__constraints = ''
        
        inp.dft = cpmd.Section()
        inp.dft.functional__blyp = '' # TO DO: update to new XC_DRIVER code
        
        # set no of steps
        inp.cpmd.maxsteps = nsteps
        
        # set timestep
        inp.cpmd.timestep = round(dt/hartree_to_ps) # convert to atomic units and round off
        
        self.inp = inp
        
        self.toYaml()
        
        return inp
