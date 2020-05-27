# This is part of MiMiCPy

"""

This module contains the handles to prepare the
MM topology, MPT and QM region/CPMD script

"""

from ..parsers import pdb as hpdb
from ..scripts import mdp
from .selector import Selector
from .base import BaseHandle
from .._global import _Global as _global
from . import _qmhelper
from ..parsers.mpt import Reader as MPTReader, write as mptwrite
from ..parsers import gro as GROReader
from ..utils.constants import hartree_to_ps, bohr_rad
from ..scripts import cpmd
from collections import defaultdict
import pandas as pd

class MM(BaseHandle):
    """
    Prepare MM topology, by running pdb2gmx, editconf, solvate, etc.
    Inherits from .core.base.BaseHandle
    
    """
    
    def __init__(self, status=[], protein=None):
        """Class constructor"""
        
        # parameters to pass to gmx genion
        self._ion_kwargs = {'pname': 'NA', 'nname': 'CL', 'neutral': ''}
        self._box_kwargs = {'c': '', 'd': 1.9, 'bt': 'cubic'} # parameters to pass to gmx editconf
        self._solavte_kwargs = {'cs': 'spc216.gro'} # parameters to pass to gmx solvate
        self._topol_kwargs = {'water': 'tip3p', 'ff': 'amber99sb-ildn'} # parameters to pass to gmx pdb2gmx
        self.his_str = '' # string version of list of histidine protonation states, input to pdb2gmx
        super().__init__(status) # call BaseHandle constructor to init _status dict
        _global.logger.write('debug', 'Set handle directory as prepareMM..')
        
        if self._status['prepMM'] == '':
            self.dir = 'prepareMM' # dir of handle, can be changed by user
        else:
            self.dir = self._status['prepMM']
        
        # all files names used when running gmx, can be changed by user
        self.confin = "confin.pdb"
        self.conf = 'conf.pdb'
        self.conf1 = 'conf1.gro'
        self.conf2 = 'conf2.gro'
        self.conf3 = 'conf3.gro'
        self.topol = "topol.top"
        self.mpt = "topol.mpt"
        self.preproc = "topol.pptop"
        self.ions = "ions"
        
        self.protein = protein
    
    def topolParams(self, **kwargs):
        """
        
        Set self._topol_kwargs, custom parameters for pdb2gmx
        Also, stringifies histiidine protonation states list
        
        """
        
        # check if his = ['D1', 'E2', ...] was passed in kwargs
        if 'his' in kwargs:
            self.his_str = '\n'.join(kwargs['his']) # join list as string
            # the list contains protonation states as human readable D1, E2
            # convert to input gromacs expects, D1-> 0; E2-> 1
            self.his_str = self.his_str.replace('D1', '0')
            self.his_str = self.his_str.replace('E2', '1')
            kwargs['his'] = '' # set the his option on, so -his will be passed to pdb2gmx
        
        self._topol_kwargs = kwargs # set custom topol kwargs
        
    def _pdb2gmx(self):
        """Runs gmx pdb2gmx for pure protien, and adds ligands/waters to pdb and .top in the right order"""
        
        self.setcurrent(key='prepMM') # set the current prepMM directory to self.dir in _status dict
        
        _global.logger.write('info', 'Preparing protein topology..')
        
        _global.logger.write('debug2', f"Writing {self.protein.name} to {self.confin}..")
        
        _global.host.write(self.protein.pdb+self.protein.water, f'{self.dir}/{self.confin}')
        
        _global.logger.write('info', 'Calculating protein topology')
        
        # call pdb2gmx, depending on if his was set
        if self.his_str == '':
            _global.logger.write('debug', "Letting Gromacs calculate histidine protonantion states..")
            self.gmx('pdb2gmx', f = self.confin, o = self.conf, dirc=self.dir, **self._topol_kwargs)
        else:
            _global.logger.write('debug', "Reading histidine protonation states from list..")
            self.gmx('pdb2gmx', f = self.confin, o = self.conf, dirc=self.dir, **self._topol_kwargs, stdin=self.his_str)
        
        conf = f'{self.dir}/{self.conf}'
        
        pdb = _global.host.read(conf) # output of pdb2gmx
        
        _global.logger.write('debug2',  f"Output of pdb2gmx saved in {self.conf}")
        
        ######
        # Followinng section reads output of pdb2gmx, which contains only protein+water residues
        # extracts crystal structure water residues
        # combine protine residues+ligand+water (IN THAT ORDER!!) are rewrite output of pdb2gmx
        ######
        lines = []
        splt = pdb.splitlines()
        for i, line in enumerate(splt[::-1]):
            vals = hpdb.readLine(line) # parse PDB using system._hndlpdb functions
            if vals['record'] != 'HETATM' and vals['record'] != 'ATOM':
                lines.append(line) 
            elif vals['record'] == 'HETATM' or vals['record'] == 'ATOM':
                if vals['resName'] == 'HOH': lines.append(line)
        
        _global.logger.write('info', "Combining ligand structure and topology with protein..")
        # combine protein, ligand and water in that order, so that gromacs doesn't complain about order
        conf_pdb = '\n'.join(splt[:len(splt)-i]) + '\n' + self.protein.ligand_pdb + '\n'.join(lines[::-1])
        
        _global.host.write(conf_pdb, conf)
        
        ######
        # Followinng section combines .itp data of all ligands into a single file
        # gromacs doesn't likes many #include for many ligands
        # so we have to combine all into a single ligands.itp file
        # Note: we need to combine all [ atomtypes ] section of itp
        # and write it as the first section of ligand.itp
        # protien.ligands is an OrderedDict so order in which we add ligands to topology
        # (i.e., order in which we loop), is also the order in which the ligands were added to pdb above
        # this is imp!! the orders have to be the same
        ######
        
        # init atomtypes section of itp
        top1 = (f"; Topology data for all non-standard resiudes in {self.protein.name} created by MiMiCPy\n"
            "; AmberTools was used to generate topolgy parameter for Amber Force Field, conversion to GMX done using Acpype"
                    "\n\n[ atomtypes ]\n")
        top2 = '' # init rest of itp
        
        _global.logger.write('debug2', "Combining ligands topology into single file ligands.itp..")
        for ligname, lig in self.protein.ligands.items():
            t1, t2 = lig.splitItp() # splits ligand itp into contents of [ atomtypes ], and eveything else
	
            top1 += t1 # update atomtypes section
            top2 += t2 # update rest of itp
            
            _global.logger.write('debug2', f"Writing position restraint file for {lig.name}")
            
            _global.host.write(lig.posre, f"{self.dir}/posre_{lig.name}.itp") # write position restraint file
            
            # add ligand to [ molecule ] section, assumed to be last section of .top
            _global.host.run(f'echo {lig.name} {lig.chains} >> topol.top')
        
        _global.host.write(top1+'\n'+top2, f"{self.dir}/ligands.itp")
        
        topol = f"{self.dir}/{self.topol}"
        
        # add #include "ligands.itp" after #include "..../forcefield.itp"
        _global.host.run(r'sed -i -r "/^#include \".+.ff\/forcefield.itp\"/a #include \"ligands.itp\"" '+topol)
        # SOL is in between protein and ligand in [ molecule ]
        # we need to remove that and add it to the end, to match pdb order
        _global.host.run('grep -v SOL topol.top > topol_.top && mv topol_.top '+topol)
        _global.host.run(f'echo SOL {self.protein.hoh_mols} >> '+topol)
        _global.logger.write('debug', "ligands.itp added to topol.top")
        
        _global.logger.write('info', 'Topology prepared..')
        
        self.saveToYaml()
        
    def ionParams(self, **kwargs): self._ion_kwargs = kwargs
    def boxParams(self, **kwargs): self._box_kwargs = kwargs
    def solvateParams(self, **kwargs): self._solavte_kwargs = kwargs
    
    def _box(self, genion_mdp):
        """Run gmx editconf, gmx solvate, gmx grompp and gmx genion in that order"""
        
        _global.logger.write('info', "Preparing system box..")
        
        self.gmx('editconf', f = self.conf, o = self.conf1, dirc=self.dir, **self._box_kwargs)
        
        _global.logger.write('info', "Solvating box..")
        self.gmx('solvate', cp = self.conf1, o = self.conf2, p = self.topol, dirc=self.dir, **self._solavte_kwargs)
        
        _global.logger.write('info', "Adding ions to neutralize charge..")
        self.grompp(genion_mdp, self.ions, gro = self.conf2, pp = self.preproc, dirc=self.dir)
        
        # sent SOL to stdin, so gromacs replaces some SOL molecules with ions 
        self.gmx('genion', s = self.ions_tpr, o = self.conf3, p = self.topol, dirc=self.dir, **self._ion_kwargs, stdin="SOL")
        
        _global.logger.write('info', 'Simulation box prepared..')
        
        self.saveToYaml()
    
    def getTopology(self, genion_mdp=mdp.MDP.defaultGenion()):
        self._pdb2gmx()
        self._box(genion_mdp)
        
        nonstd_atm_types = {}
        for ligname, lig in self.protein.ligands.items():
            nonstd_atm_types.update( dict(zip(lig.atm_types, lig.elems)) )
        
        self.getMPT(nonstd_atm_types)
    
    def getMPT(self, nonstd_atm_types={}, preproc=None, mpt=None):
        """Get the MPT topology, used in prepare.QM"""
        
        # add reading elements from protein ligands
        
        pp = self.getcurrentNone(preproc, 'pptop')
        
        if mpt == None: mpt = f"{self.dir}/{self.mpt}" # if no mpt file was passed, use default value
        
        mptwrite(pp, mpt, nonstd_atm_types)
        
        self.toYaml()
        
class QM(BaseHandle):
    """
    Prepare QM input script, by reading MPT file
    Inherits from .core.base.BaseHandle
    
    """
    
    def __init__(self, status=defaultdict(list), selector=Selector(), mpt=None, gro=None):
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
    
    def getInp(self, mdp, inp=cpmd.Input()):
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
        inp.mimic = cpmd.Section()
        inp.mimic.paths = "1\n---" #path will be set in MiMiC run function
        inp.mimic.box = '  '.join([str(s/bohr_rad) for s in self._mm_box])
        
        inp = _qmhelper.getOverlaps_Atoms(sorted_qm, inp) # get the overlaps and add atoms section of inp
        
        inp.system.charge = round(sum(self.qmatoms['charge']), 1) # system section already created in getQverlap_Atoms()
        # TO DO: give option of changing round off precision
        
        # set cpmd section
        if not inp.checkSection('cpmd'): inp.cpmd = cpmd.Section()
        inp.cpmd.mimic = ''
        inp.cpmd.parallel__constraints = ''
        
        # set no of steps from mdp file
        if not mdp.hasparam('nsteps'): mdp.nsteps = 1000 # default value, if not present in mdp
        inp.cpmd.maxsteps = mdp.nsteps
        
        # set timestep from mdp file
        if not mdp.hasparam('dt'): mdp.dt = 0.0001 # default value, if not present in mdp
        inp.cpmd.timestep = round(mdp.dt/hartree_to_ps) # convert to atomic units and round off
        
        self.inp = inp
        
        self.toYaml()
        
        return inp
