import pygmx.system.handlePDB as hpdb
from pygmx.run import mdp, gmxrun
from pygmx import host

class Prepare(gmxrun.GMX):
    
    def __init__(self, protein=None):
        print(f"Attaching protein {protein.name}..")
        self.protein = protein
        self._ion_kwargs = {'pname': 'NA', 'nname': 'CL', 'neutral': ''}
        self._box_kwargs = {'c': '', 'd': 1.9, 'bt': 'cubic'}
        self._solavte_kwargs = {'cs': 'spc216.gro'}
        
    def prepareTopology(self, histidine=[]):
        print('**Preparing protein topology**')
        
        print(f"Writing {self.protein.name} pdb to confin.pdb..")
        host.cmd.write(self.protein.pdb+self.protein.water, "confin.pdb")
        
        self.setcurrent('coords', "confin.pdb")
        
        if histidine == []:
            print("Letting Gromacs calculate histidine protonantion states..")
            self.gmx('pdb2gmx -f confin.pdb -o conf.pdb -water tip3p -ff amber99sb-ildn')
        else:
            print("Reading histidine protonation states from list..")
            his_str = '\n'.join(histidine)
            his_str = his_str.replace('D1', '0')
            his_str = his_str.replace('E2', '1')
            self.gmx('pdb2gmx', f = 'confin.pdb', o = 'conf.pdb', water = 'tip3p', ff = 'amber99sb-ildn', his='', stdin=his_str)
        
        self.setcurrent('coords', "conf.pdb")
        
        pdb = host.cmd.read('conf.pdb')
        
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
        
        host.cmd.write(conf_pdb, "conf.pdb")
        
        self.setcurrent('coords', "conf.pdb")
        
        top1 = f'; Topology data for all non-standard resiude in {self.protein.name} created by pygmx\n[ atomtypes ]\n'
        top2 = ''
        
        print("Combining lignads topology into single file ligands.itp..")
        for lig in self.protein.ligands:
            t1, t2 = lig.splitItp()
	
            top1 += t1
            top2 += t2 
            
            print(f"Writing position restraing file for {lig.name}")
            
            host.cmd.write(lig.posre, f"posre_{lig.name}.itp")
            
            host.cmd.run(f'echo {lig.name} {lig.chains} >> topol.top')
        
        host.cmd.write(top1+'\n'+top2, f"ligands.itp")
        
        host.cmd.run(r'sed -i -r "/^#include \".+.ff\/forcefield.itp\"/a #include \"ligands.itp\"" topol.top')
        host.cmd.run('grep -v SOL topol.top > topol_.top && mv topol_.top topol.top')
        host.cmd.run(f'echo SOL {self.protein.hoh_mols} >> topol.top')
        print("ligands.itp added to topol.top")
        
        self.setcurrent('topology', 'topol.top')
        
        print('**Done**')
        
    def setIons(self, **kwargs): self._ion_kwargs = kwargs
    def setBox(self, **kwargs): self._box_kwargs = kwargs
    def setSolvation(self, **kwargs): self._solavte_kwargs = kwargs
    
    def prepareBox(self, genion_mdp = mdp.MDP.defaultGenion()):
        print("**Preparing system for minimization**")
        print("Solavting protein..")
        
        self.gmx('editconf', f = self.getcurrent('coords'), o = 'conf1.gro', **self._box_kwargs)
        self.setcurrent('coords', 'conf1.gro')
        
        self.gmx('solvate', cp = 'conf1.gro',o = 'conf2.gro', p = 'topol.top', **self._solavte_kwargs)
        self.setcurrent('coords', 'conf2.gro')
        
        host.cmd.write(str(genion_mdp), "ions.mdp")
        
        print("Adding ions to neutralize charge..")
        self.gmx('grompp', f = 'ions.mdp', c = 'conf2.gro', p = 'topol.top', o = 'ions.tpr')
        
        self.gmx('genion', s = 'ions.tpr', o = 'conf3.gro', p = 'topol.top', **self._ion_kwargs, stdin="SOL")
        
        self.setcurrent('coords', 'conf3.gro')
        self.setcurrent('tpr', 'ions.tpr')
        
        print("**Done**")