# This is part of MiMiCPy

"""

This module contains helper function for prepare.QM to handle QM region prep
pptop() and getOverlaps_Atoms() adapted from the prepare_qmmm python script by Viacheslav Bolnykh

"""

from ..parsers.mpt import _mpt_writer
from ..utils.constants import bohr_rad
from ..scripts.cpmd import Atom
from collections import OrderedDict 
from ..scripts import cpmd

def _cleanqdf(qdf):
    columns = _mpt_writer.AtomsParser.columns       
    columns.extend(['x','y','z'])
    lst = [l for l in qdf.columns if l not in columns]
    return qdf.drop(lst, axis=1)

def index(qmids, name):
    """Write list of atoms to index, and returns as sting"""
    index = f'[ {name} ]\n' # name of atoms group
    for i in qmids:
        index += f'{i} '
        
    return index

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
        idx = rows['id']
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
    qm_box.extend([0, 0, 0])
    inp.system.cell = '  '.join([str(s) for s in qm_box])
    
    return inp
