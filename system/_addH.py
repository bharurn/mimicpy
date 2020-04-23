from pygmx.system import handlePDB
import pygmx.host as shell
from rdkit.Chem import PandasTools
from io import StringIO
import re

def _cleanPdb(sdf, pdb, resname):
    print("Assigning correct atom names..")
    
    atm_names = [a.GetProp('molFileAlias') for i in range(len(sdf)) for a in sdf.loc[i]['ROMol'].GetAtoms()]
    
    pdb_str = ""
    atm_i = 0 # counter for atom name
    h_i = 1 # counter for new hydrogen atom name
    # hydrogen atom names are addedas H1, H2, H3, ....
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

    return pdb_str

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
    
def do(sdf, pH, name):
    print("**Adding Hydrogens**")
    
    sdf_textio = StringIO()            
    PandasTools.WriteSDF(sdf, sdf_textio)
    
    sdf_text = sdf_textio.getvalue()
    
    print(f"Cleaning sdf..")
    sdf_text = re.sub(r'/A    1/(\n.*?)*?/M  END/', '', sdf_text, flags=re.MULTILINE)
    print("Converting to pdb using Openbabel")
    pdb = shell.cmd.run('obabel -isdf -opdb', stdin=sdf_text)
    print(f'Adding hydrogens at pH={pH}')
    pdb = shell.cmd.run(f'obabel -ipdb -p {pH} -opdb', stdin=pdb)
    
    pdb, chains = _multChains(_cleanPdb(sdf, pdb, name))
    
    pdb = shell.cmd.run('grep HETATM\|ATOM', stdin=pdb)
    print("**Done**")
    return pdb, chains

def _donoH(sdf, name):
    sdf_textio = StringIO()            
    PandasTools.WriteSDF(sdf, sdf_textio)
    
    sdf_text = sdf_textio.getvalue()
    
    print(f"Cleaning sdf..")
    sdf_text = re.sub(r'/A    1/(\n.*?)*?/M  END/', '', sdf_text, flags=re.MULTILINE)
    print("Converting to pdb using Openbabel")
    
    pdb = shell.cmd.run('obabel -isdf -opdb', stdin=sdf_text)
    
    pdb, chains = _multChains(_cleanPdb(sdf, pdb, name))
    
    pdb = shell.cmd.run('grep HETATM\|ATOM', stdin=pdb)
    print("**Done**")
    return pdb, chains