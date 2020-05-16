import mimicpy._global as _global

def _cleanprep(mol, prep_to_pdb):
    prep = _global.host.read(f'{mol}.prep')
    
    for i in prep_to_pdb:
        prep = prep.replace(i, prep_to_pdb[i])

    return prep

def _cleanItp(itp, mol):
    atoms = False
    atomtypes = False
    notwrite = False
    itp_str = ""

    mol_ = f"{mol}_" # MOL_ string
    
    hvy = []
    
    for line in itp.splitlines()[1:]:
    
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
        
        if "[ defaults ]" in line or "[ system ]" in line or "[ molecules ]" in line:
            notwrite = True
        
        if notwrite and line.strip() == '': notwrite = False
        
        if notwrite == False:
            itp_str += line+'\n' # collate line into itp_string
 
    return itp_str, hvy

def _getposre(hvy):
    posre = "[ position_restraints ]"
    posre = "[ position_restraints ]\n" + '\n'.join([f"  {i}   1  1000 1000 1000" for i in hvy])
    return posre

def do(mol, conv):
    print("Generating topology for the ligand..")
    
    prep = f"{mol}.prep"
    
    _global.host.checkFile(prep)
    
    if conv != {}:
        print(f"Mapping prep of {mol}.prep atom names to that of pdb using dictionary provided..")
        
        _global.host.write(_cleanprep(mol, conv), "params.prep")
        
        print("Ouput saved to params.prep")
        
        prep = "params.prep"
        
    print(f"Running AmberTools parmchk on {prep}..")
    
    log = _global.host.run(f'{_global.parmchk} -i {prep} -f prepi -o params.frcmod')
    
    print("Output saved to params.frcmod..")
    
    tleap_in = "source leaprc.gaff\n"
    tleap_in += f"loadamberprep {prep}\n"
    tleap_in += f"loadamberparams params.frcmod\n"
    
    if _global.host.checkFile(f"{mol}.frcmod", throw=False):
        tleap_in += f"loadamberparams {mol}.frcmod\n"
    
    tleap_in += f"saveamberparm {mol} {mol}.prmtop {mol}.inpcrd\nquit"
    
    print("Running AmberTools LEaP..")
    output = _global.host.run(f'{_global.tleap} -f -', stdin=tleap_in)
    
    log += '\n'+ tleap_in + '\n' + output + '\n'
    
    if _global.host.fileExists(f'{mol}.prmtop') == False or _global.host.fileExists(f'{mol}.inpcrd') == False:
        raise Exception(f'LEaP error!\n{output}')
    
    print(f"Output saved to {mol}.prmtop and {mol}.inpcrd..")
    
    print("Converting to Gromacs topology using Acpype..")
    output = _global.host.run(f'{_global.acpype} -p {mol}.prmtop -x {mol}.inpcrd')
    
    log += output + '\n'
    
    if _global.host.fileExists(f'{mol}_GMX.gro') == False or _global.host.fileExists(f'{mol}_GMX.top') == False:
        raise Exception(f'Acpype error!\n{output}')
        
    print(f"Output saved to {mol}_GMX.gro and {mol}_GMX.top..")
    
    print("Cleaning Acpype output..")
    itp = _global.host.read(f'{mol}_GMX.top')
    
    print("Renaming atoms to avoid conflict..")
    itp, hvy = _cleanItp(itp, mol) # rename atoms, remove uneeded sections, get list of heavy atoms for posre
    
    print("Writing position restraint data..")
    posre = _getposre(hvy)
    
    itp += f'\n; Include Position restraint file\n#ifdef POSRES\n#include "posre_{mol}.itp"\n#endif'
    
    return itp, posre, log

    