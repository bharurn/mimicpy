from ._constants import  element_names
from ..utils.scripts.cpmd import Atom
from collections import OrderedDict 
from ._constants import bohr_rad
from ..utils.scripts import cpmd

def parse_selec(selection, df):
    #selection eg., resName is SER and number < 25 and chainID not B
    ev = ''
    i = 0
    for s in selection.split():
        if i == 0:
            ev += f"(df['{s}']"
        elif s == 'or':
            ev += f' | '
            i = -1
        elif s == 'and':
            ev += f' & '
            i = -1
        elif s == 'is':
            ev += '=='
        elif s == 'not':
            ev += '!='
        elif s == '>' or s == '>=' or s == '<' or s == '<=':
            ev += s
        else:
            if s.isnumeric():
                ev += f"{s})"
            else:
                ev += f"'{s}')"
                
        i += 1

    ev = f"df.loc[{ev}]"
    ev = ev.replace("df['number']","df.index")
    return eval(ev)

def index(qmids):
    index = '[ QMatoms ]\n'
    for i in qmids:
        index += f'{i} '
        
    return index
    
def pptop(unk_lst, all_names, top):
    atm_names_to_symb = dict((el,None) for el in unk_lst)
    atm_names_to_q = dict((a,None) for a in all_names)
    
    atm_types_to_symb = {}
    start1 = False
    start2 = False
    
    for line in top.splitlines():
        if '[ atomtypes ]' in line and atm_types_to_symb == {}: # read only first atomtypes
            start1 = True
        elif start1:
            if '[' and ']' in line:
                start1 = False
            elif not line.startswith(';') and line.strip() != '':
                e = line.split()[0]
                n = int(line.split()[1])-1
                atm_types_to_symb[e]  = element_names[n]
        elif '[ atoms ]' in line:
            start2 = True
        elif start2:
            if '[' and ']' in line:
                start2 = False
            elif not line.startswith(';') and line.strip() != '':
                atm_name = line.split()[4]
                atm_type = line.split()[1]
                q = line.split()[6]
                
                if atm_name in atm_names_to_q.keys():
                    if atm_names_to_q[atm_name] == None:
                        atm_names_to_q[atm_name] = float(q)
                    
                    if atm_name in atm_names_to_symb.keys() and atm_names_to_symb[atm_name] == None:
                        atm_names_to_symb[atm_name] = atm_types_to_symb[atm_type]
        
    return atm_names_to_symb, atm_names_to_q

def getOverlaps_Atoms(qmatoms, inp):
    inp.atoms = OrderedDict()
    
    out = str(len(qmatoms))+'\n'
    mx = [None, None, None]
    mi = [None, None, None]
    for i, rows in qmatoms.iterrows():
        elem = rows['element']
        idx = rows['number']
        coords = [rows['x'], rows['y'], rows['z']]
        link = rows['link']
        
        out += f"2 {idx} 1 {i + 1}\n"
        
        if link: elem += '*'
        
        if elem not in inp.atoms.keys():
            if elem != 'H' or elem !='H*': inp.atoms[elem] = Atom(coords=[], lmax='p')
            else: inp.atoms[elem] = Atom(coords=[])
        inp.atoms[elem].coords.append([float(v)/bohr_rad for v in coords])
        
        for i, coord in enumerate(coords):
            c = float(coord)
            if mx[i] == None or c > mx[i]: mx[i] = c
            if mi[i] == None or c < mi[i]: mi[i] = c
        
    inp.mimic.overlaps = out
    inp.system = cpmd.Section()
    
    qm_box = list(map(lambda x,y: (x - y + 0.6)/bohr_rad, mx, mi))
    
    qm_box[1] = round(qm_box[1]/qm_box[0], 2)
    qm_box[2] = round(qm_box[2]/qm_box[0], 2)
    qm_box[0] = round(qm_box[0])

    inp.system.cell = '  '.join([str(s) for s in qm_box])
    return inp