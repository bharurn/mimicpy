from ..utils import handlePDB
from .._global import host as shell, obabel
import re

def _cleanPdb(sdf, pdb):
    print("Assigning correct atom names..")
    
    #####Read sdf file#####
    l_cnt = 0
    start = False
    resname = ""
    atm_names = [] # store atoms names
    elems = [] # store element list for MiMiC
    for line in sdf.splitlines():
        if l_cnt > 3 and not start: # read coordinates
            elems.append(line.split()[3]) #4th value if the element symbol
            l_cnt += 1
        elif "A    1" in line: # start reading atom names
            start = True
            val = line.split()
            if len(val) == 1:
                # get atom names to be added to the pdb file
                atm_names.append(val[0])
            elif "M  END" in line:
                start = False # stop reading atom names
            elif start:
                val = line.split()
                if len(val) == 1:
                    # get atom names to be added to the pdb file
                    atm_names.append(val[0])
                # search for ChemCompId tag to get ligand name
            elif "ChemCompId" in line:
                resname = "start"
            elif resname == "start":
                resname = line.strip()
        elif '$$$$' in line:
            l_cnt = 0 # reset l_cnt to read elements again
            

    #####END#####
    
    pdb_str = ""
    atm_i = 0 # counter for atom name
    h_i = 1 # counter for new hydrogen atom names
    # hydrogen atom names are added as H1, H2, H3, ....
    for line in pdb.splitlines():
        header = line.split()[0]
        if header == "HETATM":
            line = handlePDB.editLine(line, name=atm_names[atm_i], resName=resname)
            atm_i += 1
        elif header == "ATOM":
        #ATOM is hydrogen atom, so add H1, H2...
            hstr = "H" + str(h_i) # create hydrogen atom name
            line = handlePDB.editLine(line, name=hstr, resName=resname) # add the res name also
            h_i += 1
        
        pdb_str += line + '\n'

    return pdb_str, resname, elems

def _multChains(pdb):
    pdb_str = ""
    n_chains = 0
    stack = []
    stack_idx = 0
    
    print("Assigning correct chain IDs..")
    
    for line in pdb.splitlines():
        splt = handlePDB.readLine(line)
        
        if splt['record'] == "COMPND":
            n_chains += 1
            stack_idx = 0
            
        if splt['record'] == "HETATM" or splt['record'] == "ATOM":
            
            line = handlePDB.editLine(line, chainID=chr(ord('A')+n_chains-1), resSeq=str(n_chains))
            if splt['element'] == 'H':
                if n_chains == 1:
                    stack.append(handlePDB.readLine(line)['name'])
                else:
                    line = handlePDB.editLine(line, name=stack[stack_idx])
                    stack_idx += 1
            
        pdb_str += line + '\n'

    return pdb_str, n_chains
    
def do(sdf, pH):
    sdf_text = sdf.copy()
    
    print(f"Cleaning sdf..")
    sdf_text = re.sub(r'/A    1/(\n.*?)*?/M  END/', '', sdf_text, flags=re.MULTILINE)
    print("Converting to pdb using Openbabel")
    pdb = shell.cmd.run('{obabel} -isdf -opdb', stdin=sdf_text)
    print(f'Adding hydrogens at pH={pH}')
    pdb = shell.cmd.run(f'{obabel} -ipdb -p {pH} -opdb', stdin=pdb)
    
    pdb, resname, elems = _cleanPdb(sdf, pdb)
    pdb, chains = _multChains(pdb)
    
    pdb = shell.cmd.run('grep HETATM\|ATOM', stdin=pdb)
    
    return pdb, resname, chains, elems

def donoH(sdf):
    sdf_text = sdf.copy()
    
    print(f"Cleaning sdf..")
    sdf_text = re.sub(r'/A    1/(\n.*?)*?/M  END/', '', sdf_text, flags=re.MULTILINE)
    print("Converting to pdb using Openbabel")
    
    pdb = shell.cmd.run('{obabel} -isdf -opdb', stdin=sdf_text)
    
    pdb, resname, elems = _cleanPdb(sdf, pdb)
    pdb, chains = _multChains(pdb)
    
    pdb = shell.cmd.run('grep HETATM\|ATOM', stdin=pdb)
    return pdb, resname, chains, elems