import pyshell.local as cmd

def _cleanprep(mol, prep_to_pdb):
    file = open(f'{mol}.prep',mode='r')
 
    prep = file.read()

    for i in prep_to_pdb:
        prep = prep.replace(i, prep_to_pdb[i])

    return prep

def _renameAtoms(itp, mol):
    atoms = False
    atomtypes = False
    itp_str = ""

    mol_ = f"{mol}_" # MOL_ string
    
    hvy = []
    
    for line in itp.splitlines():
    
        if "[ atomtypes ]" in line:
            atomtypes = True # start atomtypes
        elif atomtypes:
            if line.strip() == '': # end atomtypes
                atomtypes = False
            elif not line.strip()[0] == ';': # if line is not comment
                splt = line.split() # split line
                #and stich it back together with splt[0] and splt[1], which are the atom names
                #repalced by $1_atom-names
                line = ' '+ mol_ + splt[0] + ' '*7 + mol_ + splt[1] + line[12:]
            
        if "[ atoms ]" in line:
            atoms = True # start atoms
        elif atoms:
            if line.strip() == '' or line == '[ bonds ]':
                atoms = False # end atoms
            elif not line.strip()[0] == ';' :
                splt = line.split()
                #same as above
                line = line[:9] + mol_ + splt[1] + line[11:]
                
                if 'H' not in splt[4].upper():
                    hvy.append(splt[0])
    
        itp_str += line+'\n' # collate line into itp_string
 
    return itp_str, hvy

def _getposre(hvy):
    posre = "[ position_restraints ]"
    posre = "[ position_restraints ]\n" + '\n'.join([f"  {i}   1  1000 1000 1000" for i in hvy])
    return posre

def do(mol, conv, tleap_dump=False):
    print("**Generating topology**")
    
    prep = f"{mol}.prep"
    
    cmd.checkFile(prep)
    
    if conv != {}:
        print(f"Mapping prep of {mol}.prep atom names to that of pdb..")
        
        f = open("params.prep", "w")
        f.write(_cleanprep(mol, conv))
        f.close()
    
        print("Ouput saved to params.prep")
        
        prep = "params.prep"
        
    print(f"Running AmberTools parmchk2 on {prep}..")
    
    cmd.run('parmchk2', '-i', prep, '-f', 'prepi', '-o', 'params.frcmod')
    
    print("Output saved to params.frcmod..")
    
    tleap_in = "source leaprc.gaff\n"
    tleap_in += f"loadamberprep {prep}\n"
    tleap_in += "loadamberparams params.frcmod\n"
    
    if cmd.checkFile(f"{mol}.frcmod", throw=False):
        tleap_in += f"loadamberparams {mol}.frcmod\n"
    
    tleap_in += f"saveamberparm {mol} {mol}.prmtop {mol}.inpcrd\nquit"
    
    print("Running AmberTool LEaP..")
    output = cmd.run('tleap', '-f -', stdin=tleap_in, decode=True)
    
    if tleap_dump:
        print("Dumping LEaP output:\n\n")
        print(output)
        
    print(f"Output saved to {mol}.prmtop and {mol}.inpcrd..")
    
    print("Converting to Gromacs topology using Acpype..")
    cmd.run('acpype', '-p', f'{mol}.prmtop', '-x', f'{mol}.inpcrd')
    print(f"Output saved to {mol}_GMX.gro and {mol}_GMX.top..")
    
    print("Cleaning Acpype output..")
    itp = cmd.runinSeq(['sed', '-e', '/\[ defaults \]/{N;N;N;d;}', f'{mol}_GMX.top'],\
                 ['sed', '-e', '/\[ system \]/{N;N;d;}'],\
                 ['sed', '-e', '/\[ molecules \]/{N;N;N;N;d;}'],\
                 ['sed', f's/{mol}_GMX/{mol}/g'],\
    ['sed', f's/{mol}.top created by acpype/Topology for {mol} created on/g'], decode=True)
    
    print("Renaming atoms to avoid conflict..")
    itp, hvy = _renameAtoms(itp, mol)
    
    print("Writing position restraint..")
    posre = _getposre(hvy)
    
    itp += f'\n; Include Position restraint file\n#ifdef POSRES\n#include "posre_{mol}.itp"\n#endif'
    
    print("**Done**")
    return itp, posre

    