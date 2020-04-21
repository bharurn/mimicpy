from pygmx import system
import pygmx.host as shell
from rdkit.Chem import PandasTools
import pygmx.system.handlePDB as hpdb

class Ligand:
    
    def __init__(self, pdb, itp, posre, chains, name):
        self.pdb = pdb
        self.itp = itp
        self.posre = posre
        self.chains = chains
        self.name = name
    
    def splitItp(self):
        start = False
        val = ''
        val2 = ''
        for line in self.itp.splitlines():

            if line.strip() == '[ atomtypes ]':
                start = True
            elif start:
			    # stop reading atomtypes section when blank line or moleculetype is reached
                if line.strip() == '' or  line.strip() == '[ moleculetype ]':
                    start = False
                else: val += line + '\n'
		
            elif not start: # if not start, read line into val2
                val2 += line + '\n'
		
        return val, val2	
    
    def _matchpdb2itp(self):
        start = False
        vals = {}
        i = 1
        
        print("**Matching atom order in pdb to itp**")

        print("Reading itp atoms..")
        for line in self.itp.splitlines(): 
            if "[ atoms ]" in line:
                start = True
            elif start:
                splt = line.split()
                if len(splt) == 0:
                    break
                if splt[0].isdigit():
                    vals[splt[4]] = i
                    i = i + 1
        
        pdb_str = ""
        natms = 0
        
        print("Renumbering pdb atoms..")
        for line in self.pdb.splitlines():
            splt = hpdb.readLine(line)
            if splt['record'] == "HETATM" or splt['record'] == "ATOM":
                no = vals[splt['name']]
                
                if no > natms and splt['chainID'] == 'A': natms = no
        
                no_str = str( no + natms*( ord(splt['chainID'])-ord('A') ) )
        
                line = hpdb.editLine(line, serial=no_str)
        
            pdb_str += line+'\n'
        
        self.pdb = cmd.runinSeq(['sed', "s/ATOM  /HETATM/g"],
                                ['print', "Sorting pdb atoms.."],
                                ['sort', '-nk2'], stdin=pdb_str, decode=True)
        print("**Done**")
    
    @classmethod    
    def loadFromDF(cls, mol, pH=7, prep2pdb={}, tleap_dump=False):
        mol_name = mol.loc[0]['ChemCompId']
        
        print(f"**Creating ligand {mol_name} from SDF Dataframe**")
        
        pdb, chains = system._addH.do(mol, pH, mol_name)
        
        if pH != 0:
            itp, posre = system._getItp.do(mol_name, prep2pdb, tleap_dump)
    
            lig = cls(pdb, itp, posre, chains, mol_name)
        
            lig._matchpdb2itp()
        else:
            lig = cls(pdb, "", "", 1, mol_name)
        
        return lig
    
    @classmethod    
    def loadFromFile(cls, mol_name, pH=7, prep2pdb={}, tleap_dump=False):
        cmd.checkFile(f'{mol_name}.sdf')
        mol = PandasTools.LoadSDF(f'{mol_name}.sdf')
        
        print(f"**Creating ligand {mol_name} from SDF file**")
        
        pdb, chains = system._addH.do(mol, pH, mol_name)
        
        if pH != 0:    
            itp, posre = system._getItp.do(mol_name, prep2pdb, tleap_dump)
    
            lig = cls(pdb, itp, posre, chains, mol_name)
        
            lig._matchpdb2itp()
        else:
            lig = cls(pdb, "", "", 1, mol_name)
            
        return lig
    
    @staticmethod    
    def runG09FromSDF(server, mol, nc, pH=7, parallel=False):
        
        name = mol.loc[0]['ChemCompId']
        
        print(f"**Running Gaussian Optimization for {name} on SDF Dataframe**")
        
        pdb, chains = system._addH.do(mol, pH, name)
        
        with open(f'{name}.pdb', 'w') as f: f.write(pdb)
        
        print("Generating Gaussian input file using AmberTools Antechamber..")
        
        cmd.run("antechamber", "-i", f"{name}.pdb", "-fi", "pdb", "-o", f"{name}.com", "-fo", "gcrt", "-nc", f"{nc}")
        
        print(f'Transferring files to {server.name}..')
        
        server.sftp.put(f'{name}.com', f'{name}.com', confirm=True)
        
        print(f'Running Gaussian on {server.name}..')
        
        server.query_rate = 30
        server.redirectStdout(f'{name}.out')
        
        out = server.run(f'g09 {name}.com {name}.out')
        
        print("**Done**")
        
        return out
    
    @staticmethod
    def getPrepfromG09(mol):
        print(f"**Converting Gaussian output file to Amber prep**")
        cmd.checkFile(f"{mol}.out")
        print('Running AmberTools antechamber..')
        out = cmd.run("antechamber", "-i", f"{mol}.out", "-fi", "gout", "-o", f"{mol}.prep", "-fo", "prepi", "-c", "resp", "-rn", mol)
        print(f"Antechamber output dumped...\nResidue {mol} created in {mol}.prep..\n**Done**")
        return out
        
class StdResidue(Ligand):
    @classmethod
    def loadFromDF(cls, mol):
        return super().loadFromDF(mol, 0)
    
    @classmethod 
    def loadFromFile(cls, mol):
        return super().loadFromFile(mol, 0)
    