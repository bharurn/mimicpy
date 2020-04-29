from . import _addH, _getItp
from ..utils import handlePDB as hpdb
from .. import host as shell
from rdkit.Chem import PandasTools

class Ligand:
    
    def __init__(self, pdb, itp, posre, chains, name, elems):
        self.pdb = pdb
        self.itp = itp
        self.posre = posre
        self.chains = chains
        self.name = name
        self.elems = elems
    
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
        
        self.pdb = pdb_str.replace('ATOM  ','HETATM')
        print("Sorting pdb atoms..")
        self.pdb = shell.cmd.run('sort -nk2', stdin=self.pdb)
        print("**Done**")
    
    @classmethod
    def _load(cls, mol, pH, mol_name, prep2pdb, tleap_dump):
        pdb, chains = _addH.do(mol, pH, mol_name)
        
        itp, posre = _getItp.do(mol_name, prep2pdb, tleap_dump)
        
        elems = [a.GetSymbol() for a in mol.loc[0]['ROMol'].GetAtoms()] # for MiMiC
        
        lig = cls(pdb, itp, posre, chains, mol_name, elems)
        
        lig._matchpdb2itp()
        
        return lig
    
    @classmethod    
    def loadFromDF(cls, mol, pH=7, prep2pdb={}, tleap_dump=False):
        mol_name = mol.loc[0]['ChemCompId']
        
        print(f"**Creating ligand {mol_name} from SDF Dataframe**")
        
        return cls._load(mol, pH, mol_name, prep2pdb, tleap_dump)
    
    @classmethod    
    def loadFromFile(cls, mol_name, pH=7, prep2pdb={}, tleap_dump=False):
        shell.cmd.checkFile(f'{mol_name}.sdf')
        
        mol = PandasTools.LoadSDF(f'{mol_name}.sdf')
        
        print(f"**Creating ligand {mol_name} from SDF file**")
        
        return cls._load(mol, pH, mol_name, prep2pdb, tleap_dump)
    
    @staticmethod    
    def runG09FromSDF(server, mol, nc, pH=7, parallel=False):
        
        name = mol.loc[0]['ChemCompId']
        
        print(f"**Running Gaussian Optimization for {name} on SDF Dataframe**")
        
        pdb, chains = _addH.do(mol, pH, name)
        
        with open(f'{name}.pdb', 'w') as f: f.write(pdb)
        
        print("Generating Gaussian input file using AmberTools Antechamber..")
        
        shell.cmd.run(f"antechamber -i {name}.pdb -fi pdb -o {name}.com -fo gcrt -nc {nc}")
        
        print(f'Running Gaussian on {server.name}..')
        
        shell.cmd.query_rate = 30
        shell.cmd.redirectStdout(f'{name}.out') # add this to local shell
        
        out = shell.cmd.run(f'g09 {name}.com {name}.out')
        
        print("**Done**")
        
        return out
    
    @staticmethod
    def getPrepfromG09(mol):
        print(f"**Converting Gaussian output file to Amber prep**")
        shell.cmd.checkFile(f"{mol}.out")
        print('Running AmberTools antechamber..')
        out = shell.cmd.run(f"antechamber -i {mol}.out -fi gout -o {mol}.prep -fo prepi -c resp -rn {mol}")
        print(f"Antechamber output dumped...\nResidue {mol} created in {mol}.prep..\n**Done**")
        
        return out
        
class StdResidue(Ligand):
    @classmethod
    def loadFromDF(cls, mol):
        mol_name = mol.loc[0]['ChemCompId']
        
        print(f"**Creating standard residue {mol_name} from SDF Dataframe*")
        
        pdb, chains = _addH._donoH(mol, mol_name)
        return cls(pdb, "", "", 1, mol_name)
    
    @classmethod 
    def loadFromFile(cls, mol_name):
        shell.cmd.checkFile(f'{mol_name}.sdf')
        mol = PandasTools.LoadSDF(f'{mol_name}.sdf')
        
        print(f"**Creating standard residue {mol_name} from SDF file**")
        
        pdb, chains = _addH._donoH(mol, mol_name)
        
        return cls(pdb, "", "", 1, mol_name)
    