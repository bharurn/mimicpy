from . import ligand as lig, _hndlWater
from . import _hndlpdb as hpdb
import mimicpy._global as _global
import requests
from collections import defaultdict, OrderedDict

class Protein:
   
   def __init__(self, pdb, name):
        self.name = name
        
        print(f"Creating protein {name}..")
        
        self.pdb = pdb
        
        print("Extracting water..")
        
        self.water = _global.host.run('grep ^HETATM', stdin=self.pdb)
        self.water = _global.host.run('grep HOH', stdin=self.water)
        
        print("Extracting amino acid residues..")
        
        self.pdb  = f"TITLE     {self.name}\n" + _global.host.run('grep ^ATOM', stdin=self.pdb)
        
        self.ligands=OrderedDict()
        
        self.ligand_pdb = ""
        
        self._lig_elems = [] # for MiMiC
        
   def addLigand(self, ligand):
        
       if not isinstance(ligand, lig.NonStdLigand):
           raise Exception(f"Cannot add {type(ligand)} as ligand")
       
       if not isinstance(ligand, lig.StdLigand):
           print(f"Adding non standard ligand {ligand.name} to {self.name}..") 
           
           self.ligands.update({ligand.name:ligand})
           
           self.ligand_pdb += ligand.pdb
           
           self._lig_elems.extend(ligand.elems) # for MiMiC
           
       else:
           print(f"Adding standard ligand {ligand.name} to {self.name}..") 
           
           self.pdb += ligand.pdb
    
   def addLigands(self, ligand_list):
       for ligand in ligand_list:  self.addLigand(ligand)
       
   def stripWater(self, *args):
       
      ids = [ r for (res, chain, dist) in args\
                             for r in _hndlWater.coordsWithin(self.ligands[res].pdb, chain, self.water, dist) ]
     
      self.hoh_mols = len(ids)
      
      command = '\|'.join(ids)
    
      self.water = _global.host.run(f'grep {command}', stdin=self.water)
    
   @classmethod
   def fromRCSB(cls, pdbid, chains=None, howToreturn=0):
       print(f"Accessing PDB {pdbid} from RCSB database..")
       print("Downloading PDB file..")
           
       response = requests.get(f'http://files.rcsb.org/view/{pdbid}.pdb')
       response.raise_for_status()
       
       pdb = response.text
       
       ligs = defaultdict(str)
       
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

        
                ligs[v[0]] += resp.text
                
            elif vals['record'] == 'HETNAM': break
            
       pdb_ = ''
       
       for line in pdb.splitlines():
            vals = hpdb.readLine(line)
        
            if vals['record'] == 'ATOM' or vals['record'] == 'HETATM':
                if vals['chainID'] in chains:
                    pdb_ += line + '\n'
       
       if howToreturn == 0:
           prt = cls(pdb_, pdbid)
           prt.addLigands(ligs.values())
           print("Returned as Protein with Ligands added..")
           return prt
       else:
           print("Returned as raw data..")
           return pdb_, ligs
       
   @classmethod
   def fromFile(cls, pdb):
        _global.host.checkFile(pdb)
        f = _global.host.read(pdb)
        return cls(f, pdb.split('.')[0])