# This is part of MiMiCPy

"""

This module contains helper function for prepare.QM to handle QM region prep
pptop() and getOverlaps_Atoms() adapted from the prepare_qmmm python script by Viacheslav Bolnykh

"""

from ..parsers import mpt
from ..utils.constants import bohr_rad, hartree_to_ps
from ..scripts.cpmd import Atom
from collections import OrderedDict 
from ..scripts import cpmd

def cleanqdf(qdf):
    columns = mpt.columns.copy() # copy it otherwise original gets edited    
    columns.extend(['x','y','z'])
    col_to_drop = [l for l in qdf.columns if l not in columns]
    
    qdf.index = qdf.index.set_names(['id']) # rename index, so it is id when we reset index in prepare.getInp()
    
    return qdf.drop(col_to_drop, axis=1)

def index(qmids, name, space_len=6, col_len=15):
    """Write list of atoms to index, and returns as sting"""
    
    max_len = len(str(max(qmids)))
    spaces = space_len if max_len <= space_len else max_len
    
    index = f'[ {name} ]\n' # name of atoms group
    for i, idx in enumerate(qmids):
        if i%col_len == 0: index += '\n'
        index += "{:{}}".format(idx, spaces)
    
    return index

def check_mdp(mdp, nsteps = 1000, dt = 0.0001):
    
    if mdp is None:
        return nsteps, dt/hartree_to_ps, ""
    
    mdp_errors = []
        
    if not mdp.hasparam('integrator') or mdp.integrator != 'mimic':
        mdp_errors.append("\tintegrator = mimic not set")
    
    if not mdp.hasparam('qmmm_grps') or mdp.qmmm_grps != 'QMatoms':
        mdp_errors.append("\tQMMM-grps = QMatoms not set")
        
    if mdp.hasparam('constraints') and mdp.constraints != 'none':
        mdp_errors.append("\tMolecule should not be constrained by Gromacs, set constraints = none")
    
    if mdp.hasparam('tcoupl') and mdp.tcoupl != 'no':
        mdp_errors.append("\tTemperature coupling will not be active, set tcoupl = no")
    
    if mdp.hasparam('pcoupl') and mdp.pcoupl != 'no':
        mdp_errors.append("\tPressure coupling will not be active, set pcoupl = no")
    
    #check for other errors in mdp file
    
    if mdp.hasparam('nsteps'):    
        nsteps = int(mdp.nsteps)
    
    if mdp.hasparam('dt'):    
        dt = float(mdp.dt)
    
    return nsteps, dt/hartree_to_ps, "\n".join(mdp_errors)

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
        
        if link: elem += '*' # append astericks to link atoms, so they are added as normal elem dict value
        
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
    if not inp.checkSection('system'): inp.system = cpmd.Section() # init system section of script
    inp.system.poisson__solver__tuckerman = ''
    inp.system.symmetry = 0

    # box size fro mx and mi
    # add 0.7 nm for Poisson solver's requirement
    qm_box = list(map(lambda x,y: (x - y + 0.7)/bohr_rad, mx, mi))
    qm_box[1] = round(qm_box[1]/qm_box[0], 1)
    qm_box[2] = round(qm_box[2]/qm_box[0], 1)
    qm_box[0] = round(qm_box[0])
    qm_box.extend([0, 0, 0])
    inp.system.cell = '  '.join([str(s) for s in qm_box])
    inp.system.cutoff = 70.0
    
    return inp
