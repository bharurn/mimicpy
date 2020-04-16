import pygmx.system.ligand as lig
import pyshell.local as cmd
from pygmx.system import _hndlWater
import requests
import pygmx.system.handlePDB as hpdb
from collections import defaultdict
import io
from rdkit.Chem import PandasTools

class Protein:
   
   def __init__(self, pdb, name):
        self.name = name
        
        print(f"**Creating protein {name}**")
        
        self.pdb = pdb
        
        print("Extracting water..")
        
        self.water = cmd.runinSeq(['grep', '^HETATM'],['grep', "HOH"], stdin=self.pdb, decode=True)
        
        print("Extracting amino acid residues..")
        
        self.pdb  = f"TITLE     {self.name}\n" + cmd.run('grep', "^ATOM", stdin=self.pdb, decode=True)
        
        self.ligands=[]
        self.ligands_dict={}
        
        self.ligand_pdb = ""
        
        print("**Done**")
        
   def addLigand(self, ligand):
        
       if not isinstance(ligand, lig.Ligand):
           raise Exception(f"Cannot add {type(ligand)} as ligand")
       
       if not isinstance(ligand, lig.StdResidue):
           print(f"Adding non standard ligand {ligand.name} to {self.name}..") 
           
           self.ligands.append(ligand)
           
           self.ligands_dict.update({ligand.name:ligand})
           
           self.ligand_pdb += ligand.pdb
       else:
           print(f"Adding standard ligand {ligand.name} to {self.name}..") 
           
           self.pdb += ligand.pdb
       
   def stripWater(self, *args):
       
      ids = [ r for (res, chain, dist) in args\
                             for r in _hndlWater.coordsWithin(self.ligands_dict[res].pdb, chain, self.water, dist) ]
     
      self.hoh_mols = len(ids)
      
      command = '\|'.join(ids)
    
      self.water = cmd.run('grep', command, stdin=self.water, decode=True)
     
   def prepare(self, ssh, sftp, gmx_load, gmx="gmx_mpi"):
       f = sftp.open('confin.pdb', 'w')
       f.write(self.pdb+'END')
       f.close()
       
       stdin,stdout,stderr=ssh.exec_command(f"module load {gmx_load}", get_pty=True)
       stdin,stdout,stderr=ssh.exec_command(f"{gmx} pdb2gmx -f confin.pdb -o conf.pdb -water tip3p -ff amber99sb-ildn", get_pty=True)
       
       
       return stdout.read(), stderr.read()
   
   @classmethod
   def loadFromRCSB(cls, pdbid, chains=None):
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
       
       return cls(pdb_, pdbid), ligs
       
   @classmethod
   def loadFromFile(cls, pdb):
        cmd.checkFile(f"{pdb}.pdb")
        
        f = open(f"{pdb}.pdb", "r")
        
        return cls(f.read(), pdb)