from ..utils import handlePDB as hpdb
from ..utils.scripts import mdp
from .base import Run
import mimicpy._global as _global
from . import _qmhelper
from ._constants import bohr_rad
from ..utils.scripts import cpmd

class MM(Run):
    
    def __init__(self, protein=None):
        print(f"Attaching protein {protein.name}..")
        self.protein = protein
        self._ion_kwargs = {'pname': 'NA', 'nname': 'CL', 'neutral': ''}
        self._box_kwargs = {'c': '', 'd': 1.9, 'bt': 'cubic'}
        self._solavte_kwargs = {'cs': 'spc216.gro'}
        super().__init__()
        
    def prepareTopology(self, histidine=[]):
        print('Preparing protein topology..')
        
        print(f"Writing {self.protein.name} pdb to confin.pdb..")
        _global.host.write(self.protein.pdb+self.protein.water, "confin.pdb")
        
        self.setcurrent('coords', "confin.pdb")
        
        if histidine == []:
            print("Letting Gromacs calculate histidine protonantion states..")
            Run.gmx('pdb2gmx -f confin.pdb -o conf.pdb -water tip3p -ff amber99sb-ildn')
        else:
            print("Reading histidine protonation states from list..")
            his_str = '\n'.join(histidine)
            his_str = his_str.replace('D1', '0')
            his_str = his_str.replace('E2', '1')
            Run.gmx('pdb2gmx', f = 'confin.pdb', o = 'conf.pdb', water = 'tip3p', ff = 'amber99sb-ildn', his='', stdin=his_str)
        
        self.setcurrent('coords', "conf.pdb")
        
        pdb = _global.host.read('conf.pdb')
        
        lines = []
        splt = pdb.splitlines()
        for i, line in enumerate(splt[::-1]):
            vals = hpdb.readLine(line)
            if vals['record'] != 'HETATM' and vals['record'] != 'ATOM':
                lines.append(line)  
            elif vals['resName'] == 'HOH':
                lines.append(line)
            else:
                break
        
        print("Adding non-standard residues to structure..")
        conf_pdb = '\n'.join(splt[:len(splt)-i]) + '\n' + self.protein.ligand_pdb + '\n'.join(lines[::-1])
        
        _global.host.write(conf_pdb, "conf.pdb")
        
        self.setcurrent('coords', "conf.pdb")
        
        top1 = (f"; Topology data for all non-standard resiudes in {self.protein.name} created by MiMiCPy\n"
            "; AmberTools was used to generate topolgy parameter for Amber Force Field, conversion done using Acpype"
                    "\n\n[ atomtypes ]\n")
        top2 = ''
        
        print("Combining ligands topology into single file ligands.itp..")
        for ligname, lig in self.protein.ligands.items():
            t1, t2 = lig.splitItp()
	
            top1 += t1
            top2 += t2 
            
            print(f"Writing position restraint file for {lig.name}")
            
            _global.host.write(lig.posre, f"posre_{lig.name}.itp")
            
            _global.host.run(f'echo {lig.name} {lig.chains} >> topol.top')
        
        _global.host.write(top1+'\n'+top2, f"ligands.itp")
        
        _global.host.run(r'sed -i -r "/^#include \".+.ff\/forcefield.itp\"/a #include \"ligands.itp\"" topol.top')
        _global.host.run('grep -v SOL topol.top > topol_.top && mv topol_.top topol.top')
        _global.host.run(f'echo SOL {self.protein.hoh_mols} >> topol.top')
        print("ligands.itp added to topol.top")
        
        self.setcurrent('topology', 'topol.top')
        print('Topology prepared..')
        
        self.saveToYaml()
        
    def setIons(self, **kwargs): self._ion_kwargs = kwargs
    def setBox(self, **kwargs): self._box_kwargs = kwargs
    def setSolvation(self, **kwargs): self._solavte_kwargs = kwargs
    
    def prepareBox(self, genion_mdp = mdp.MDP.defaultGenion()):
        print("Preparing system box..")
        print("Solavting protein..")
        
        Run.gmx('editconf', f = self.getcurrent('coords'), o = 'conf1.gro', **self._box_kwargs)
        self.setcurrent('coords', 'conf1.gro')
        
        Run.gmx('solvate', cp = 'conf1.gro',o = 'conf2.gro', p = 'topol.top', **self._solavte_kwargs)
        self.setcurrent('coords', 'conf2.gro')
        
        _global.host.write(str(genion_mdp), "ions.mdp")
        
        print("Adding ions to neutralize charge..")
        Run.gmx('grompp', f = 'ions.mdp', c = 'conf2.gro', p = 'topol.top', o = 'ions.tpr')
        
        Run.gmx('genion', s = 'ions.tpr', o = 'conf3.gro', p = 'topol.top', **self._ion_kwargs, stdin="SOL")
        
        self.setcurrent('coords', 'conf3.gro')
        self.setcurrent('tpr', 'ions.tpr')
        print('Simulation box prepared..')
        
        self.saveToYaml()
        
class QM(Run):
    
    def __init__(self, prepare):
        self._qmatoms = {}
        self.continueFrom(prepare)
        self.protein = prepare.protein
        self._protein_res = []
        self._mm_box = [5, 5, 5]
        self.inp = cpmd.Input()
        
    def addAtoms(self, **kwargs):
        
        lst = ['serial', 'name', 'resName', 'chainID', 'element']
        
        for k in kwargs.keys():
            if k not in lst:
                slst = ", ".join(lst)
                raise Exception(f'{k} not a valid selection. Please use the following keywords:\n{slst}')
        
        pdb_coords = _global.host.read(self.gethistory('coords')[1]) #pdb file after pdb2gmx, with ligands included
        gro_coords = _global.host.read(self.getcurrent('coords'))
        
        pdb_splt = pdb_coords.splitlines()
        gro_splt = gro_coords.splitlines()
        
        for i, pdb, gro in enumerate(zip(pdb_splt, gro_splt)):
            if hpdb.matchLine(pdb, kwargs):
                idx = i+1
                
                resname = hpdb.readLine(pdb)['resName'] # check for non standard ligands
                if resname in self.protein.ligands: elem = self.protein._lig_elems[i]
                else:
                    elem = self.protein._prt_atom_names[i]
                    self._protein_res.append(idx)
                
                coords = [float(v)/bohr_rad for v in gro.split()[2:5]]
                
                self._qmatoms[idx] = (elem, coords)
        
        self._mm_box = [float(v)/bohr_rad for v in gro_splt[-1].split()]
    
    def addLinkAtoms(self, **kwargs):
        pass
    
    def detectLinkAtoms(self):
        pass
    
    def prepareQMRegion(self, mdp, inp=cpmd.Input()):
        mdp.integrator = 'mimic'
        mdp.nsteps = 10000 # dummy value
        mdp.dt = 0.002 #dummy value
        mdp.QMMM_grps = 'QMatoms'
        
        _global.host.write(str(mdp), 'mimic.mdp')
        _global.host.write(_qmhelper.index(self.qmatoms.keys()), 'index.ndx')
        
        Run.gmx('grompp', f='mimic.mdp', c=self.getcurrent('coords'),\
                 p=self.getcurrent('topology'), o='mimic.tpr', pp='processed.top', n='index.ndx')
        self.setcurrent('tpr', 'mimic.tpr')
        
        qmatoms_updated = _qmhelper.pptop(self._qmatoms, self._protein_res, _global.host.read('processed.top'))
        
        inp.mimic = cpmd.Section()
        inp.mimic.paths = f"1\n{_global.host.pwd()}/gmx"
        inp.mimic.box = '  '.join(self._mm_box)
        
        inp = _qmhelper.getOverlaps_Atoms(qmatoms_updated, inp)
        
        if not inp.checkSection('cpmd'): inp.cpmd = cpmd.Section()
        inp.cpmd.mimic = ''
        inp.cpmd.parallel__constraints = ''
        
        self.inp = inp
        
        return inp