import os
import logging
import pandas as pd
from ..topology.mpt import Mpt
from ..scripts.mdp import Mdp
from ..scripts.ndx import Ndx
from ..scripts.cpmd import CpmdScript, Pseudopotential
from ..utils.errors import SelectionError
from ..utils.constants import BOHR_RADIUS
from ..utils.file_handler import write


class Preparation:

    def __init__(self, selector):
        self.__qm_atoms = pd.DataFrame()
        self.selector = selector

    @staticmethod
    def __clean_qdf(qdf):
        columns = Mpt.columns.copy()
        columns.extend(['x', 'y', 'z'])
        columns_to_drop = [l for l in qdf.columns if l not in columns]
        qdf.index = qdf.index.set_names(['id'])
        return qdf.drop(columns_to_drop, axis=1)

    def add(self, selection=None, is_link=False):
        qdf = Preparation.__clean_qdf(self.selector.select(selection))
        qdf.insert(2, 'is_link', [int(is_link)]*len(qdf))
        self.__qm_atoms = self.__qm_atoms.append(qdf)

    def delete(self, selection=None):
        qdf = Preparation.__clean_qdf(self.selector.select(selection))
        self.__qm_atoms = self.__qm_atoms.drop(qdf.index, errors='ignore')

    def clear(self):
        self.__qm_atoms = pd.DataFrame()

    @property
    def qm_atoms(self):
        return self.__qm_atoms

    def get_mimic_input(self, inp_tmp=None, mdp_inp=None, ndx_out=None, inp_out=None):
        """Args:
            inp_tmp: cpmd input file, used as template
            mdp_inp: gromacs input file, checked for errors
            ndx_out: gromacs index file, output
            inp_out: mimic cpmd input file, output
        """

        def qm_cell():
            dims = [0, 0, 0]
            for i, r in enumerate(['x', 'y', 'z']):
                dims[i] = (abs(max(self.__qm_atoms[r]) - min(self.__qm_atoms[r])) + 0.7)/BOHR_RADIUS
            a, b, c = dims
            cell = ' '.join((str(round(a, 1)), str(round(b/a, 1)), str(round(c/a, 1)), '0 0 0'))
            return cell

        # Check for obvious errors in selection
        if self.__qm_atoms.empty:
            raise SelectionError('No atoms have been selected for the QM partition.')

        # Delete self.mpt for better garbage collection
        #try:
        #    del self.selector.mpt
        #except AttributeError:
        #    pass

        # Retrieve number of steps and timestep from mdp_inp and do some checks
        if mdp_inp is None:
            maxsteps, timestep, mdp_errors = 1000, 5.0, ['Using default values for maxstep and timestep.']
        else:
            maxsteps, timestep, mdp_errors = Mdp.from_file(mdp_inp).check()
        for error in mdp_errors:
            logging.warning('%s', error)

        # Create an index group in GROMACS format (and write it to a file)
        qm_ndx_group = Ndx('qmatoms') # use default name
        qm_ndx_group.qmatoms = self.__qm_atoms.index.to_list()
        if ndx_out:
            write(str(qm_ndx_group), ndx_out, 'w')
            logging.info('Wrote Gromacs index file to %s', ndx_out)

        # Create CPMD input script
        sorted_qm_atoms = self.__qm_atoms.sort_values(by=['is_link', 'element']).reset_index()

        if inp_tmp is None:
            cpmd = CpmdScript('Cpmd', 'System', 'Mimic', 'Atoms')
        else:
            cpmd = CpmdScript.from_file(inp_tmp)

        # Get overlaps and atoms
        overlaps = f'{len(sorted_qm_atoms)}'
        for i, atom in sorted_qm_atoms.iterrows():
            gromacs_id = atom['id']
            cpmd_id = i + 1
            overlaps += f'\n2 {gromacs_id} 1 {cpmd_id}'
            element = str(atom['element']).lower()
            coords = [atom['x']/BOHR_RADIUS, atom['y']/BOHR_RADIUS, atom['z']/BOHR_RADIUS]
            if atom['is_link']:
                element += '_link'
            if cpmd.atoms.has_parameter(element):
                pp_block = getattr(cpmd.atoms, element)
                pp_block.coords.append(coords)
            else:
                setattr(cpmd.atoms, element, Pseudopotential(coords))

        if not cpmd.mimic.has_parameter('paths'):
            cpmd.mimic.paths = '1\n' + str(os.getcwd())

        cpmd.mimic.overlaps = overlaps
        cpmd.mimic.box = ' '.join([str(s/BOHR_RADIUS) for s in self.selector.mm_box])
        cpmd.system.cell = qm_cell()

        total_charge = sum(self.__qm_atoms['charge'])
        if not round(total_charge, 2).is_integer():
            logging.warning('Total charge of QM region is %s. Rounding to integer.', total_charge)
        cpmd.system.charge = round(total_charge)

        cpmd.cpmd.maxsteps = maxsteps
        cpmd.cpmd.timestep = timestep

        if inp_out is None:
            logging.info('Created new CPMD input script for MiMiC run')
        else:
            write(str(cpmd), inp_out, 'w')
            logging.info('Wrote new CPMD input script to %s', inp_out)

        return qm_ndx_group, cpmd
