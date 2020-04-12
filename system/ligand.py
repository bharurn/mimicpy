from pygmx import system
import pyshell.local as cmd
from rdkit.Chem import PandasTools
import pygmx.system.handlePDB as hpdb

class Ligand:
    
    def __init__(self, pdb, itp, posre, chains, name):
        self.pdb = pdb
        self.itp = itp
        self.posre = posre
        self.chains = chains
        self.name = name
    
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
    def load(cls, mol, pH=7, prep2pdb={}, tleap_dump=False):
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
