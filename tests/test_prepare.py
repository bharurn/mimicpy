import io
import logging
from mimicpy import Mdp
from mimicpy.utils.errors import SelectionError
import pandas as pd
import pytest

class MockSelector():
    def __init__(self):
        dct = {'type': [], 'resid': [], 'resname': [],\
            'name': [], 'charge': [], 'element': [], 'mass': [], 'mol': [], 'x': [], 'y': [], 'z': []}

        with open('prepare_test_files/qmatoms.txt', 'r') as f:
            for i, line in enumerate(f.readlines()):
                if i == 0:
                    self.box = line.split()
                elif i == 1:
                    mol = line.strip()
                elif i > 1:
                    name, x, y, z, q = line.split()
                    dct['type'].append(name)
                    dct['name'].append(name)
                    dct['element'].append(name[0])
                    dct['resname'].append(mol)
                    dct['mol'].append(mol)
                    dct['resid'].append(1)
                    dct['mass'].append(0)
                    dct['x'].append(float(x))
                    dct['y'].append(float(y))
                    dct['z'].append(float(z))
                    dct['charge'].append(float(q))

        self.df = pd.DataFrame(dct)
        self.df.index = self.df.index + 1

        self.mm_box = [float(i) for i in self.box]

    def select(self, selection=None):
         return self.df

def read_mock_cpmd():
    lst = {}
    with open('prepare_test_files/cpmd_params.txt') as f:
        lines = f.readlines()
        lst['mm_box'] = [float(i) for i in lines[0].split()]
        lst['qm_box'] = [float(i) for i in lines[1].split()]
        lst['charge'] = float(lines[2])

    with open('prepare_test_files/overlaps.txt') as f:
        lst['overlaps'] = set(f.read().splitlines()[1:])

    from collections import defaultdict
    atoms = defaultdict(list)

    with open('prepare_test_files/cpmd_atoms.xyz', 'r') as f:
        for i, line in enumerate(f.readlines()):
            if i > 1:
                name, x, y, z = line.split()
                atoms[name].append([float(x),float(y),float(z)])

    return lst, atoms

def is_close(a, b, rel_tol=1e-09, abs_tol=0.0):
    return abs(a-b) <= max(rel_tol * max(abs(a), abs(b)), abs_tol)

def is_close_array(a, b, rel_tol=1e-09, abs_tol=0.0):
    for i,j in zip(a, b):
        if not isinstance(i, float) and not isinstance(j, float):
            is_close_array(i, j, rel_tol, abs_tol)

        elif not is_close(i, j, rel_tol, abs_tol):
            return False

    return True

def test_prepare():
    warns = io.StringIO()

    fileh = logging.StreamHandler(warns)
    formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
    fileh.setFormatter(formatter)

    log = logging.getLogger()  # root logger
    for hdlr in log.handlers[:]:  # remove all old handlers
        log.removeHandler(hdlr)
    log.addHandler(fileh)      # set the new handler

    from mimicpy.core.prepare import Preparation
    prep = Preparation(MockSelector())

    with pytest.raises(SelectionError) as e:
        assert prep.get_mimic_input()
    assert str(e.value) == "No atoms have been selected for the QM partition."

    prep.clear()

    mock_mdp = Mdp(name='test', tcoupl='yes')

    prep.add('resid is 1')
    ndx, cpmd = prep.get_mimic_input(mdp_inp=mock_mdp)

    import os
    assert cpmd.mimic.paths == f"1\n{os.getcwd()}"

    lst, atoms = read_mock_cpmd()

    box = [float(i) for i in cpmd.mimic.box.split()]
    cell = [float(i) for i in cpmd.system.cell.split()]

    assert is_close_array(box, lst['mm_box'], 1e-03)
    assert is_close_array(cell, lst['qm_box'], 1e-01)
    assert cpmd.system.charge == lst['charge']

    overlaps = set(cpmd.mimic.overlaps.splitlines()[1:])
    assert overlaps == lst['overlaps']

    assert is_close_array(cpmd.atoms.c.coords, atoms['C'], 1e-03)
    assert is_close_array(cpmd.atoms.h.coords, atoms['H'], 1e-03)
    assert is_close_array(cpmd.atoms.o.coords, atoms['O'], 1e-03)

    assert "Total charge of QM region" in warns.getvalue()\
        and "Rounding to integer" in warns.getvalue()

    assert "Wrong integrator for MiMiC run, set integrator = mimic" in warns.getvalue()
    assert "Index group for QM atoms is not qmatoms, set QMMM-grps to the appropriate group" in warns.getvalue()
    assert "Temperature coupling will not be active, set tcoupl = no" in warns.getvalue()
    assert "Molecules should not be constrained by Gromacs, set constraints = none" not in warns.getvalue()
    assert "Pressure coupling will not be active, set pcoupl = no" not in warns.getvalue()