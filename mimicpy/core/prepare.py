# This is part of MiMiCPy

"""

This module contains the handles to prepare the
MM topology, MPT and QM region/CPMD script

"""

from ..system import _hndlpdb as hpdb
from ..scripts import mdp
from .base import BaseHandle
from .simulate import MD
from .._global import _Global as _global
from . import _qmhelper, _mpt_helper
from ._constants import hartree_to_ps
from ..scripts import cpmd
from collections import defaultdict

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
        
        self.dir = 'prepareMM' # dir of handle, can be changed by user
        
        # all files names used when running gmx, can be changed by user
        self.confin = "confin.pdb"
        self.conf = 'conf.pdb'
        self.conf1 = 'conf1.gro'
        self.conf2 = 'conf2.gro'
        self.conf3 = 'conf3.gro'
        self.topol = "topol.top"
        self.mpt = "topol.mpt"
        self.ions_mdp = "ions.mdp"
        self.ions_tpr = "ions.tpr"
        
        if protein:
            # if protien was passed, automatically do everything
            self.getTopology(protein)
            self.getBox()
    
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
        
    def getTopology(self, protein):
        """Runs gmx pdb2gmx for pure protien, and adds ligands/waters to pdb and .top in the right order"""
        
        self.setcurrent(key='prepMM') # set the current prepMM directory to self.dir in _status dict
        
        _global.logger.write('info', 'Preparing protein topology..')
        
        _global.logger.write('debug2', f"Writing {protein.name} to {self.confin}..")
        
        _global.host.write(protein.pdb+protein.water, f'{self.dir}/{self.confin}')
        
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
        # writes only protein residues to MPT row by row
        # read ligand pdb and adds that to MPT
        # combine protine residues+ligand+water (IN THAT ORDER!!) are rewrite output of pdb2gmx
        ######
        lines = []
        mpt = _mpt_helper.MPTWriter()
        splt = pdb.splitlines()
        for i, line in enumerate(splt[::-1]):
            vals = hpdb.readLine(line) # parse PDB using system._hndlpdb functions
            if vals['record'] != 'HETATM' and vals['record'] != 'ATOM':
                lines.append(line) 
            elif vals['record'] == 'HETATM' or vals['record'] == 'ATOM':
                if vals['resName'] == 'HOH': lines.append(line)
                else:
                    #create topology dataframe
                    mpt.write_row(vals)
        
        # read ligand pdb data and put into mpt
        for l in protein.ligand_pdb.splitlines():
            #TO DO: edit resno of ligand to match rest of pdb
            # not imp for gromacs but imp when saving as mpt
            mpt.write_row(hpdb.readLine(line))
        
        mpt.write(self.mpt, self.dir)

        _global.logger.write('info', "Combining ligand structure and topology with protein..")
        # combine protein, ligand and water in that order, so that gromacs doesn't complain about order
        conf_pdb = '\n'.join(splt[:len(splt)-i]) + '\n' + protein.ligand_pdb + '\n'.join(lines[::-1])
        
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
        top1 = (f"; Topology data for all non-standard resiudes in {protein.name} created by MiMiCPy\n"
            "; AmberTools was used to generate topolgy parameter for Amber Force Field, conversion to GMX done using Acpype"
                    "\n\n[ atomtypes ]\n")
        top2 = '' # init rest of itp
        
        _global.logger.write('debug2', "Combining ligands topology into single file ligands.itp..")
        for ligname, lig in protein.ligands.items():
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
        # SOL is in between protien and ligand in [ molecule ]
        # we need to remove that and add it to the end, to match pdb order
        _global.host.run('grep -v SOL topol.top > topol_.top && mv topol_.top '+topol)
        _global.host.run(f'echo SOL {protein.hoh_mols} >> '+topol)
        _global.logger.write('debug', "ligands.itp added to topol.top")
        
        _global.logger.write('info', 'Topology prepared..')
        
        self.saveToYaml()
        
    def ionParams(self, **kwargs): self._ion_kwargs = kwargs
    def boxParams(self, **kwargs): self._box_kwargs = kwargs
    def solvateParams(self, **kwargs): self._solavte_kwargs = kwargs
    
    def getBox(self, genion_mdp = mdp.MDP.defaultGenion()):
        """Run gmx editconf, gmx solvate, gmx grompp and gmx genion in that order"""
        
        _global.logger.write('info', "Preparing system box..")
        
        self.gmx('editconf', f = self.conf, o = self.conf1, dirc=self.dir, **self._box_kwargs)
        
        _global.logger.write('info', "Solvating box..")
        self.gmx('solvate', cp = self.conf1, o = self.conf2, p = self.topol, dirc=self.dir, **self._solavte_kwargs)
        
        _global.host.write(str(genion_mdp), f'{self.dir}/{self.ions_mdp}')
        
        _global.logger.write('info', "Adding ions to neutralize charge..")
        self.gmx('grompp', f = self.ions_mdp, c = self.conf2, p = self.topol, o = self.ions_tpr, dirc=self.dir)
        
        # sent SOL to stdin, so gromacs replaces some SOL molecules with ions 
        self.gmx('genion', s = self.ions_tpr, o = self.conf3, p = self.topol, dirc=self.dir, **self._ion_kwargs, stdin="SOL")
        
        _global.logger.write('info', 'Simulation box prepared..')
        
        self.saveToYaml()
    
    def getMPT(self, pdb=None, mpt=None):
        """Get the MPT topology, used in prepare.QM for systems where getTopology() was not run"""
        
        # Related to Issue 1: thsi reads whole PDB including water in one shot
        # very slow, change!!
        
        # assume pdb has protein and ligand info
        pdb = self.getcurrentNone(pdb, 'pdb')
        pdb_data = _global.host.read(pdb)
        
        mpt = _mpt_helper.MPTWriter()
        splt = pdb_data.splitlines()
        for line in splt:
            vals = hpdb.readLine(line)
            if vals['record'] == 'HETATM' or vals['record'] == 'ATOM':
                if vals['resName'] == 'HOH' or vals['resName'] == 'SOL': break
                else: mpt.write_row(vals)
        
        if mpt == None: mpt = self.mpt # if no mpt file was passed, use default value
        
        mpt.write(self.mpt, self.dir)
        
class QM(MD):
    """
    Prepare QM input script, by reading MPT file
    Allows for usage of human readble selection language to add atoms to QM section
    
    """
    
    def __init__(self, status=defaultdict(list), mpt=None, coords=None):
        """Class constructor"""
        
        super().__init__(status) # call BaseHandle.__init__() to init status dict
        
        # combine mpt with coordinate file, and get dataframe, also returns the last line of gro
        # TO DO: check if latest run is trr or gro, and if trr convert
        self.df, self._mm_box = _mpt_helper.read(self.getcurrentNone(mpt, 'mpt'), self.getcurrentNone(coords, 'gro'))
        
        # init scripts and paths/files
        self.inp = cpmd.Input()
        self.qmatoms = None
        self.dir = 'prepareQM'
        self.index = 'index.ndx'
        self.preprc = 'processed.top'
        self.mdp_tpr = 'mimic'
        
    def add(self, selection, link=False):
        """Add to QM region using selection langauage"""
        # call the parse selection function
        qdf = _qmhelper.parse_selec(selection, self.df)
            
        qdf.insert(2, 'link', [int(link)]*len(qdf)) # add a new column link which is integer of link argument
        
        # add qdf to self.qmatoms, append if already exists
        if self.qmatoms is None:
            self.qmatoms = qdf
        else:
            self.qmatoms = self.qmatoms.append(qdf)
    
    def delete(self, selection):
        """Delete from QM region using selection langauage"""
        remove = _qmhelper.parse_selec(selection, self.df)
        # drop selection, ignore errors to ignore extra residues selected that are not present in self.qmatoms
        self.qmatoms = self.qmatoms.drop(remove.index, errors='ignore')
    
    def clear(self):
        """Empty the QM region to start over"""
        self.qmatoms = None
    
    def getInp(self, mdp, inp=cpmd.Input()):
        """
        Create the QM region from the atoms added
        Steps done:
            Set integrator to mimic, and QMMM_grps to QMatoms in mdp script
            Writes index files of atoms in QM region with name QMatoms
            Run gmx grompp to get mimic.tpr and processed.top
            Read processed.top to fill in missing element symbols and all charge info
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
        self.grompp(mdp, f'{self.mdp_tpr}', pp=self.preprc, n=self.index, dirc=dirc)
        
        unk_lst = self.qmatoms[self.qmatoms['element'].apply(lambda x: False if x.strip() != '' else True)]['name'].to_list()
        
        _global.logger.write('debug', "Reading force field data to fill in charges and missing atomic symbol information..")
        # the full processed.top file is loaded into memory for fast manipulation, may not be feasible for large files
        conv, q_conv = _qmhelper.pptop(unk_lst, self.qmatoms['name'].to_list(), \
                                       _global.host.read(f'{dirc}/{self.preprc}') ) # read element conv and charge conv
        # element conv: dict of atom names -> element symbols, q cond: dict of atom names -> charge
       
        # iterate through df to combine name & element, get charges for each name also
        # pass lambda to df.apply to run lambda for each row
        # lambda returns a tuple of (elem_symbol, charge) for each row
        # which apply returns as a list of tuples, we then use zip() to get it as two sepearte lists
        # logic of lambda: if element column is empty: then return conv[atom_name] for element name, else just return element name
        elems, q = zip(*self.qmatoms[['name', 'element']].apply(\
                       lambda x: (x[1], q_conv[x[0]]) if x[1].strip() != '' else (conv[x[0]], q_conv[x[0]]), axis=1))
        
        # we cannot edit self.qmatoms directly, as pandas complains
        # so do it in round out way, we only assign elements, as only sum(q) is needed
        qm_df = self.qmatoms.assign(element=elems)
        self.qmatoms = qm_df # reassign it to original df
        
        # sort by link column first, then element symbol
        # ensures that all link atoms are last, and all elements are bunched together
        # index is also reset, for getOverlaps_Atoms()
        self.qmatoms = self.qmatoms.sort_values(by=['link', 'element']).reset_index()
        
        _global.logger.write('info', "Creating CPMD input script..")
        inp.mimic = cpmd.Section()
        inp.mimic.paths = "1\n---" #path will be set in MiMiC run function
        inp.mimic.box = '  '.join([str(s) for s in self._mm_box])
        
        inp = _qmhelper.getOverlaps_Atoms(self.qmatoms, inp) # get the overlaps and add atoms section of inp
        
        inp.system.charge = round(sum(q), 1) # system section already created in getQverlap_Atoms()
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