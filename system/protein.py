from . import ligand as lig, _hndlWater
from ..utils import handlePDB as hpdb
from .. import host as shell
import requests
from collections import defaultdict
import io
from rdkit.Chem import PandasTools
from collections import OrderedDict

class Protein:
   
   def __init__(self, pdb, name):
        self.name = name
        
        print(f"**Creating protein {name}**")
        
        self.pdb = pdb
        
        print("Extracting water..")
        
        self.water = shell.cmd.run('grep ^HETATM', stdin=self.pdb)
        self.water = shell.cmd.run('grep HOH', stdin=self.water)
        
        print("Extracting amino acid residues..")
        
        self.pdb  = f"TITLE     {self.name}\n" + shell.cmd.run('grep ^ATOM', stdin=self.pdb)
        
        self.ligands=OrderedDict()
        
        self.ligand_pdb = ""
        
        print("**Done**")
        
   def addLigand(self, ligand):
        
       if not isinstance(ligand, lig.Ligand):
           raise Exception(f"Cannot add {type(ligand)} as ligand")
       
       if not isinstance(ligand, lig.StdResidue):
           print(f"Adding non standard ligand {ligand.name} to {self.name}..") 
           
           self.ligands.update({ligand.name:ligand})
           
           self.ligand_pdb += ligand.pdb
       else:
           print(f"Adding standard ligand {ligand.name} to {self.name}..") 
           
           self.pdb += ligand.pdb
       
   def stripWater(self, *args):
       
      ids = [ r for (res, chain, dist) in args\
                             for r in _hndlWater.coordsWithin(self.ligands[res].pdb, chain, self.water, dist) ]
     
      self.hoh_mols = len(ids)
      
      command = '\|'.join(ids)
    
      self.water = shell.cmd.run(f'grep {command}', stdin=self.water)
    
   @classmethod
   def loadFromRCSB(cls, pdbid, chains=None, howToreturn=0):
       print(f"**Accessing PDB {pdbid} from RCSB database**")
       print("Downloading PDB file..")
           
       response = requests.get(f'http://files.rcsb.org/view/{pdbid}.pdb')
       response.raise_for_status()
       
       pdb = response.text
       
       lig = {}

       lig = defaultdict(str)
       
       print("Parsing ligands..")
       
       for l in pdb.splitlines():
            vals = hpdb.readLine(l)
            
            if vals['record'].strip() == 'HET':
                
                if chains is not None and vals['content'].split()[1] not in chains:
                    continue
                
                v = vals['content'].split()[:3]
        
                query = '_'.join(v)
                
                print(f"Downloading ligands {v[0]}, chain {v[1]}")
        
                resp = requests.get(f"https://files.rcsb.org/cci/download/{pdbid}_{query}_NO_H.sdf")
                resp.raise_for_status()

        
                lig[v[0]] += resp.text
                
            elif vals['record'] == 'HETNAM': break
    
       ligs = {}
       
       for lg in lig:
            output = io.BytesIO(lig[lg].encode('utf-8'))

            ligs[lg] = PandasTools.LoadSDF(output)
            
       pdb_ = ''
       
       for line in pdb.splitlines():
            vals = hpdb.readLine(line)
        
            if vals['record'] == 'ATOM' or vals['record'] == 'HETATM':
                if vals['chainID'] in chains:
                    pdb_ += line + '\n'
            
       print("**Done**")
       
       if howToreturn == 0:
           return cls(pdb_, pdbid), ligs
       else:
           return pdb_, ligs
       
   @classmethod
   def loadFromFile(cls, pdb):
        shell.cmd.checkFile(f"{pdb}.pdb")
        
        f = shell.cmd.read(f"{pdb}.pdb")
        
        return cls(f, pdb)