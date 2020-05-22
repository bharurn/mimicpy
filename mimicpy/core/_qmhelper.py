# This is part of MiMiCPy

"""

This module contains helper function for prepare.QM to handle QM region prep
pptop() and getOverlaps_Atoms() adapted from the prepare_qmmm python script by Viacheslav Bolnykh

"""

from ._constants import  element_names
from ..scripts.cpmd import Atom
from collections import OrderedDict 
from ._constants import bohr_rad
from ..scripts import cpmd
from .._global import _Global as gbl

def parse_selec(selection, df):
    """Translate selection language string into pandas dataframe selection"""
    
    # if selection is a lambda func, then just call and return it
    # this is provided for debuggin puposes
    LAMBDA = lambda:0
    if isinstance(selection, type(LAMBDA)):
        return df[selection(df)]
            
    # below code translates selection langauge to pandas boolean
    # selection eg., resName is SER and number < 25 and chainID not B
    # will be translated to df['resName'] == 'SER' and df.index == 25 and df['chainID'] != 'B'
    
    ev = '' # cint onverted string
    i = 0 # counter to keep track of word position
    for s in selection.split():
        if i == 0: # if starting of set
            ev += f"(df['{s}']"
        # if and/or encountered, reset i to -1 (will become 0 due to i+= 1 at end)
        # so we can start parsing again
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
        else: # parse everything else, meant for the third word
            if s.isnumeric():
                ev += f"{s})"
            else:
                ev += f"'{s}')"
                
        i += 1

    ev = f"df.loc[{ev}]" # eg., df.loc[ df['resName'] == 'SER' and df.index == 25 and df['chainID'] != 'B' ]
    ev = ev.replace("df['number']","df.index") # replace df['number'] to df.index as number is the index of df
    gbl.logger.write('debug2', f'Selection command translated to: {ev}')
    return eval(ev) # evaluate string and return the dataframe

def index(qmids, name):
    """Write list of atoms to index, and returns as sting"""
    index = f'[ {name} ]\n' # name of atoms group
    for i in qmids:
        index += f'{i} '
        
    return index
    
def pptop(unk_lst, all_names, top):
    """
    Parse prerocessed topology and return mapping of 
    atom_name->element symbol & atom_name->parital charge
    
    Logic:
        Read first [ atomtypes ] section and get mapping from atom types->element symbol
        All other [ atomtypes ] are not read as these correspond to non std ligands and we already have the symbol info
        for those from system.ligands.NonStdLigands
        Read all [ atoms ] sections:
            Get mapping of atom names->q
            Get mapping of atom names->element symbol through atom types->element symbol mapping
    """
    
    atm_names_to_symb = dict((el,None) for el in unk_lst) # init atm_name->symb dict with only those that are unknown
    atm_names_to_q = dict((a,None) for a in all_names) # init atm_name->q dict
    
    atm_types_to_symb = {} # init atom types to symbol
    start1 = False
    start2 = False
    
    for line in top.splitlines():
        if '[ atomtypes ]' in line and atm_types_to_symb == {}: # read only first atomtypes
            start1 = True
        elif start1:
            if '[' and ']' in line:
                start1 = False
            elif not line.startswith(';') and line.strip() != '': # ignore comments and null lines
                e = line.split()[0] # first val is atom type
                n = int(line.split()[1])-1 # second is element no.
                atm_types_to_symb[e]  = element_names[n] # fill up atom types->symbol mapping
        elif '[ atoms ]' in line: # read all atom sections
            start2 = True
        elif start2:
            if '[' and ']' in line:
                start2 = False
            elif not line.startswith(';') and line.strip() != '': # ignore comments and null lines
                atm_name = line.split()[4] # fifth val is atom name
                atm_type = line.split()[1] # second val is atom type
                q = line.split()[6] # seventh val is charge
                
                if atm_name in atm_names_to_q.keys():
                    if atm_names_to_q[atm_name] == None:
                        atm_names_to_q[atm_name] = float(q) # atom names to charge
                    
                    if atm_name in atm_names_to_symb.keys() and atm_names_to_symb[atm_name] == None:
                        atm_names_to_symb[atm_name] = atm_types_to_symb[atm_type] # atoms names to symbol
        
    return atm_names_to_symb, atm_names_to_q

def getOverlaps_Atoms(qmatoms, inp):
    """
    Fill up ATOMS section of CPMD script with qmatoms
    The dataframe is assumed to be ordered correctly
    The qm boz size is also calculated
    atom_name->element symbol & atom_name->parital charge
    
    Logic:
        Read first [ atomtypes ] section and get mapping from atom types->element symbol
        All other [ atomtypes ] are not read as these correspond to non std ligands and we already have the symbol info
        for those from system.ligands.NonStdLigands
        Read all [ atoms ] sections:
            Get mapping of atom names->q
            Get mapping of atom names->element symbol through atom types->element symbol mapping
    """
    
    inp.atoms = OrderedDict() # init atoms as OrderedDict, since order is important
    # the keys are element symbols, and values are scripts.cpmd.Atoms() objects
    
    out = str(len(qmatoms))+'\n' # init overlap section string with no of atoms
    mx = [None, None, None] # max coords
    mi = [None, None, None] # min coords
    for i, rows in qmatoms.iterrows():
        elem = rows['element']
        idx = rows['number']
        coords = [rows['x'], rows['y'], rows['z']]
        link = rows['link']
        
        out += f"2 {idx} 1 {i + 1}\n" # overlap section string
        
        if link: elem += '*'
        
        if elem not in inp.atoms.keys(): # if new element is found
            # init a new Atom() for that element key
            if elem != 'H' or elem !='H*': inp.atoms[elem] = Atom(coords=[], lmax='p')
            else: inp.atoms[elem] = Atom(coords=[], lmax='s')
        # append coords to coords variable of Atom() object of that element
        inp.atoms[elem].coords.append([float(v)/bohr_rad for v in coords])
        
        for i, coord in enumerate(coords): # find max and min coords in all 3 directions
            c = float(coord)
            if mx[i] == None or c > mx[i]: mx[i] = c
            if mi[i] == None or c < mi[i]: mi[i] = c
        
    inp.mimic.overlaps = out
    inp.system = cpmd.Section() # init system section of script
    
    # box size fro mx and mi
    qm_box = list(map(lambda x,y: (x - y + 0.6)/bohr_rad, mx, mi))
    qm_box[1] = round(qm_box[1]/qm_box[0], 2)
    qm_box[2] = round(qm_box[2]/qm_box[0], 2)
    qm_box[0] = round(qm_box[0])
    inp.system.cell = '  '.join([str(s) for s in qm_box])
    
    return inp