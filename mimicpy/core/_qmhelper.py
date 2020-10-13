
from ..io import mpt
from ..utils.constants import bohr_rad, hartree_to_ps
from ..scripts.cpmd import Atom
from collections import OrderedDict
from ..scripts import cpmd

def cleanqdf(qdf):
    columns = mpt.columns.copy()
    columns.extend(['x','y','z'])
    columns_to_drop = [l for l in qdf.columns if l not in columns]
    qdf.index = qdf.index.set_names(['id'])
    return qdf.drop(columns_to_drop, axis=1)

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

