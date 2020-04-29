from . import _addH, _getItp
from ..utils import handlePDB as hpdb
from ..utils._misc import f
from .._global import host as shell, gauss

class NonStdLigand: 
    
    def __init__(self, pdb, itp, posre, chains, name, elems):
        self.pdb = pdb
        self.itp = itp
        self.posre = posre
        self.chains = chains
        self.name = name
        self.elems = elems # for MiMiC
    
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
        
        print("Matching atom order in pdb to itp..")

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
    
    @classmethod
    def _load(cls, mol, pH, prep2pdb, tleap_dump):
        pdb, mol_name, chains, elems = _addH.do(mol, pH)
        
        itp, posre = _getItp.do(mol_name, prep2pdb, tleap_dump)
        
        lig = cls(pdb, itp, posre, chains, mol_name, elems)
        
        lig._matchpdb2itp()
        
        return lig
    
    @classmethod    
    def loadFromDF(cls, mol, pH=7, prep2pdb={}, tleap_dump=False):
        print(f"Creating non standard ligand from SDF data..")
        
        return cls._load(mol, pH, prep2pdb, tleap_dump)
    
    @classmethod    
    def loadFromFile(cls, file, pH=7, prep2pdb={}, tleap_dump=False):
        shell.checkFile(f'{file}.sdf')
        
        mol = shell.read(file)
        
        print(f"Creating non standard ligand from {file}..")
        
        return cls._load(mol, pH, prep2pdb, tleap_dump)
    
    @staticmethod    
    def runGaussFromSDF(mol, nc, pH=7, parallel=False):
        
        print(f"Running Gaussian Calculation on ligand..")
        
        print(f"Converting to PDB..")
        
        pdb, chains, name, elems = _addH.do(mol, pH)
        
        shell.write(pdb, f(name,'.pdb'))
        
        print("Generating Gaussian input file using AmberTools Antechamber..")
        
        shell.run(f("antechamber -i ",name,".pdb -fi pdb -o ",name,".com -fo gcrt -nc ",nc))
        
        shell.query_rate = 30
        shell.redirectStdout(f(name,'.out')) # add this to local shell
        
        if gauss is None:
            raise Exception("No Gaussian executable given!")
        
        out = shell.cmd.run(f(gauss,' ',name,'.com ',name,'.out'))
        
        return out
    
    @staticmethod
    def getPrepfromGauss(mol):
        print(f"Converting Gaussian output file to Amber prep..")
        shell.cmd.checkFile(f(mol,".out"))
        print('Running AmberTools antechamber..')
        out = shell.cmd.run(f"antechamber -i {mol}.out -fi gout -o {mol}.prep -fo prepi -c resp -rn {mol}")
        print(f"Antechamber output dumped...\nResidue {mol} created in {mol}.prep..\n**Done**")
        
        return out
        
class StdLigand(NonStdLigand):
    @classmethod
    def loadFromDF(cls, mol):
        print(f"Creating standard ligand from SDF data..")
        
        pdb, chains, mol_name, elems = _addH.donoH(mol)
        return cls(pdb, "", "", 1, mol_name, mol_name, elems)
    
    @classmethod 
    def loadFromFile(cls, file):
        shell.checkFile(f'{file}.sdf')
        
        mol = shell.read(file)
        
        print(f"Creating standard ligand from {file}..")
        
        pdb, chains, mol_name, elems = _addH.donoH(mol)
        return cls(pdb, "", "", 1, mol_name, mol_name, elems)
    