from ._constants import  element_names
from ..utils.scripts.cpmd import Atom
from collections import OrderedDict

def index(qmids):
    index = '[ QMatoms ]\n'
    for i in qmids:
        index += f'{i} '
        
    return index
    
def pptop(qm, prt_res, top):
    atm_names_to_symb = {q[0]:None for i,q in enumerate(qm) if i in prt_res}
    atm_types_to_symb = {}
    start1 = False
    start2 = False
    
    for line in top.splitlines():
        if '[ atomtypes ]' in line:
            start1 = True
        elif start1:
            if not line.startswith(';') and line.strip() != '':
                e = line.split()[0]
                n = line.split()[1]
                atm_types_to_symb[e] = element_names[n]
            elif '[' and ']' in line:
                start1 = False
        elif '[ atoms ]' in line:
            start2 = True
        elif start2:
            if not line.startswith(';') and line.strip() != '':
                atm_name = line.split()[4]
                atm_type = line.split()[1]
                
                if atm_name in atm_names_to_symb.keys():
                    if atm_names_to_symb[atm_name] == None:
                        atm_names_to_symb[atm_name] = atm_types_to_symb[atm_type]
    
    qm_order_by_symb = {}
                    
    for i in prt_res:
        atm_name = qm[i][0]
        qm_order_by_symb[atm_names_to_symb[atm_name]] = (i, qm[i][1])
        
    return OrderedDict(sorted(qm_order_by_symb.items()))

def getOverlaps_Atoms(qmatoms, inp):
    for e in set(qmatoms.keys()):
        inp.atoms[e] = Atom() # initialize atom sections with empty atoms
    
    out = str(len(qmatoms))
    for i, elem in enumerate(qmatoms.keys()):
        idx = qmatoms[elem][0]
        out += f"2 {idx} 1 {i + 1}\n"
        coords = qmatoms[elem][1]
        inp.atoms[elem].coord.append(coords)
        
    inp.mimic.overlaps = out
    return inp