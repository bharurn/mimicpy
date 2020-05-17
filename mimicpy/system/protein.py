from . import ligand as lig, _hndlWater
from . import _hndlpdb as hpdb
from .._global import _Global as _global
import urllib.request as req
from collections import defaultdict, OrderedDict

class Protein:
   
   def __init__(self, pdb, name):
        self.name = name
        
        _global.logger.write('debug', f"Creating protein {name}..")
        
        self.pdb = pdb
        
        _global.logger.write('debug', "Extracting water..")
        
        self.water = _global.host.run('grep ^HETATM', stdin=self.pdb)
        self.water = _global.host.run('grep HOH', stdin=self.water)
        
        _global.logger.write('debug', "Extracting amino acid residues..")
        
        self.pdb  = f"TITLE     {self.name}\n" + _global.host.run('grep ^ATOM', stdin=self.pdb)
        
        self.ligands=OrderedDict()
        
        self.ligand_pdb = ""
        
        self._lig_elems = [] # for MiMiC
        
   def addLigand(self, ligand):
        
       if not isinstance(ligand, lig.NonStdLigand):
           raise TypeError(f"Cannot add {type(ligand)} as ligand")
       
       if not isinstance(ligand, lig.StdLigand):
           _global.logger.write('debug', f"Adding non standard ligand {ligand.name} to {self.name}..") 
           
           self.ligands.update({ligand.name:ligand})
           
           self.ligand_pdb += ligand.pdb
           
           self._lig_elems.extend(ligand.elems) # for MiMiC
           
       else:
           _global.logger.write('debug', f"Adding standard ligand {ligand.name} to {self.name}..") 
           
           self.pdb += ligand.pdb
    
   def addLigands(self, ligand_list):
       for ligand in ligand_list:  self.addLigand(ligand)
       
   def stripWater(self, *args):
       # args: list of ['lig_name', 'chain', dist] eg. ['AKG', 'A', 3], ['AKG', 'B', 3]
      ids = [ r for (res, chain, dist) in args\
                             for r in _hndlWater.coordsWithin(self.ligands[res].pdb, chain, self.water, dist) ]
     
      self.hoh_mols = len(ids)
      
      command = '\|'.join(ids)
    
      self.water = _global.host.run(f'grep {command}', stdin=self.water)
    
   @classmethod
   def fromRCSB(cls, pdbid, chains=None, howToreturn=0):
       _global.logger.write('info', f"Accessing PDB {pdbid} from RCSB database..")
       _global.logger.write('debug', "Downloading PDB file..")
          
       url = f'http://files.rcsb.org/view/{pdbid}.pdb'
       _global.logger.write('debug', f"Openeing {url}..")
       response = req.urlopen(url)
       #response.raise_for_status()
       
       pdb = response.read().decode('utf-8')
       
       ligs = defaultdict(str)
       
       _global.logger.write('info', "Downloading ligands..")
       
       for l in pdb.splitlines():
            vals = hpdb.readLine(l)
            
            if vals['record'].strip() == 'HET':
                
                if chains is not None and vals['content'].split()[1] not in chains:
                    continue
                
                v = vals['content'].split()[:3]
        
                query = '_'.join(v)
                
                _global.logger.write('debug', f"Downloading ligands {v[0]}, chain {v[1]}")
        
                url = f"https://files.rcsb.org/cci/download/{pdbid}_{query}_NO_H.sdf"
                _global.logger.write('debug', f"Openeing {url}..")
                resp = req.urlopen(url)
                #resp.raise_for_status()

                ligs[v[0]] += resp.read().decode('utf-8')
                
            elif vals['record'] == 'HETNAM': break
            
       pdb_ = ''
       
       for line in pdb.splitlines():
            vals = hpdb.readLine(line)
            
            if vals['record'] == 'ATOM' or vals['record'] == 'HETATM':
                if chains is None or vals['chainID'] in chains:
                    pdb_ += line + '\n'
       
       if howToreturn == 0:
           prt = cls(pdb_, pdbid)
           for v in ligs.values():
               prt.addLigands(lig.NonStdLigand(v))
           _global.logger.write('info', "Returned as protein object with ligands added..")
           return prt
       else:
           _global.logger.write('info', "Returned as raw data..")
           return pdb_, ligs
       
   @classmethod
   def fromFile(cls, pdb):
        _global.host.checkFile(pdb)
        f = _global.host.read(pdb)
        return cls(f, pdb.split('.')[0])