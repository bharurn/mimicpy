import pygmx.system.handlePDB as hpdb
from pygmx.run import mdp, gmxrun

class Prepper(gmxrun.GMXRun):
    
    def __init__(self, server, protein, modules=[], sources=[], gmx="gmx", directory=".", ignloadxerr = True):
        print("**Setting up Prepper**")
        print(f"Attaching protein {protein.name}..")
        self.protein = protein
        super().__init__(server, modules, sources, gmx, directory, ignloadxerr)
        
        self.ic = ''
        self.pion = 'NA'
        self.nion = 'CL'
        
    def prepareTopology(self, histidine=[]):
        print('**Preparing system topology**')
        
        print(f"Writing {self.protein.name} pdb to confin.pdb..")
        with self.shell.vi("confin.pdb", 'w') as f:
            f.write(self.protein.pdb+self.protein.water)
        
        if histidine == []:
            print("Letting Gromacs calculate histidine protonantion states..")
            self.shell.run(f'gmx_mpi pdb2gmx -f confin.pdb -o conf.pdb -water tip3p -ff amber99sb-ildn')
        else:
            print("Reading histidine protonation states from list..")
            his_str = '\n'.join(histidine)
            his_str = his_str.replace('D1', '0')
            his_str = his_str.replace('E2', '1')
            self.shell.run(f'gmx_mpi pdb2gmx -f confin.pdb -o conf.pdb -water tip3p -ff amber99sb-ildn -his',\
                           stdin=his_str)
        
        print("Reading conf.pdb..")
        with self.shell.vi("conf.pdb", 'r') as f:
            pdb = f.read().decode('utf-8')
        
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
        
        print("Adding non-standard residues to conf.pdb..")
        conf_pdb = '\n'.join(splt[:len(splt)-i]) + '\n' + self.protein.ligand_pdb + '\n'.join(lines[::-1])
        
        with self.shell.vi("conf.pdb", 'w') as f:
            f.write(conf_pdb)
        
        top1 = f'; Topology data for all non-standard resiude in {self.protein.name} created by pygmx\n[ atomtypes ]\n'
        top2 = ''
        
        print("Combining lignads topology into single file ligands.itp..")
        for lig in self.protein.ligands:
            t1, t2 = lig.splitItp()
	
            top1 += t1
            top2 += t2 
            
            print(f"Writing position restraing file for {lig.name}")
            f = self.shell.vi(f"posre_{lig.name}.itp", 'w')
            f.write(lig.posre)
            f.close()
        
            self.shell.run(f'echo {lig.name} {lig.chains} >> topol.top')
        
        with self.shell.vi(f"ligands.itp", 'w') as f:
            f.write(top1+'\n'+top2)
        
        self.shell.run(r'sed -i -r "/^#include \".+.ff\/forcefield.itp\"/a #include \"ligands.itp\"" topol.top')
        self.shell.run('grep -v SOL topol.top > topol_.top && mv topol_.top topol.top')
        self.shell.run(f'echo SOL {self.protein.hoh_mols} >> topol.top')
        print("ligands.itp added to topol.top")
        
        print('**Done**')
        
    def setIons(self, pion=None, nion=None, ion_conc=None):
        if ion_conc: self.ic = f"-conc {ion_conc}"
        if pion: self.pion = pion
        if nion: self.nion = nion
    
    def prepareForEM(self, genion_mdp = mdp.MDP.defaultGenion(), em_mdp = mdp.MDP.defaultEM()):
        print("**Preparing system for minimization**")
        print("Reading conf.pdb and outputting conf1.gro..")
        self.shell.run('gmx_mpi editconf -f conf.pdb -o conf1.gro -c -d 1.9 -bt cubic')
        
        self.shell.run('gmx_mpi solvate -cp conf1.gro -cs spc216.gro -o conf2.gro -p topol.top')
        print("Output saved to conf2.gro..")
        
        with self.shell.vi("ions.mdp", 'w') as f:
            f.write(str(genion_mdp))
        
        print("Adding ions to neutralize charge..")
        self.shell.run('gmx_mpi grompp -f ions.mdp -c conf2.gro -p topol.top -o ions.tpr')
        
        self.shell.run(f'gmx_mpi genion -s ions.tpr -o conf3.gro -p topol.top -pname {self.pion} -nname {self.nion} -neutral {self.ic}',\
                       stdin="SOL")
        
        print("Output saved to conf3.gro..")
        
        with self.shell.vi("em.mdp", 'w') as f:
            f.write(str(em_mdp))
        
        print("Generating portable run file for minimization..")
        
        self.shell.run('gmx_mpi grompp -f em.mdp -c conf3.gro -p topol.top -o em.tpr')
        
        print("Writting output to prepper.log..")
        
        with open('prepper.log', 'w')  as f: f.write(self.shell.log)
        
        print("**Done**")