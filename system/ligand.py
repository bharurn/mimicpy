from pygmx import system
import pyshell.local as cmd
from rdkit.Chem import PandasTools

class Ligand:
    
    def __init__(self, pdb, itp, posre, chains, name):
        self.pdb = pdb
        self.itp = itp
        self.posre = posre
        self.chains = chains
        self.name = name
    
    @classmethod    
    def loadFromDF(cls, mol, pH=7, prep2pdb={}, tleap_dump=False):
        mol_name = mol.loc[0]['ChemCompId']
        
        print(f"**Creating ligand {mol_name} from SDF Dataframe**")
        
        pdb, chains = system._addH.do(mol, pH, mol_name)
        
        if pH != 0:
            itp, posre = system._getItp.do(mol_name, prep2pdb, tleap_dump)
    
            lig = cls(pdb, itp, posre, chains, mol_name)
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
