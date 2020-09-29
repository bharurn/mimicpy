import logging
import pandas as pd
from .selector import GroSelector
from ..io.mpt import Mpt
from ..io.gro import Gro
from ..scripts.mdp import Mdp
from ..scripts.cpmd import InputScript, Section
from .._global import _Global as gbl
from ..utils.errors import SelectionError
from ..utils.constants import BOHR_RADIUS


class Preparation:

    def __init__(self, topology, structure):
        self.mpt = Mpt.from_file(topology)  # TODO: Rename from_file since topology can also be a mpt object
        self.structure = Gro(structure)  # TODO: Write adapter for several formats
        self.qm_atoms = pd.DataFrame()
        self.selector = GroSelector(self.mpt, self.structure)  # TODO: Write adapter for other selectors

    @staticmethod
    def __clean_qdf(qdf):
        columns = Mpt.columns.copy()
        columns.extend(['x', 'y', 'z'])
        columns_to_drop = [l for l in qdf.columns if l not in columns]
        qdf.index = qdf.index.set_names(['id'])
        return qdf.drop(columns_to_drop, axis=1)

    @staticmethod
    def __ndx_group(qm_ids, group_name):
        col_len = 15
        space_len = 6
        max_len = len(str(max(qm_ids)))
        spaces = space_len if max_len <= space_len else max_len
        index = f'[ {group_name} ]'
        for i, idx in enumerate(qm_ids):
            if i%col_len == 0:
                index += '\n'
            index += "{:{}}".format(idx, spaces) # TODO: Stick to f'...' syntax
        return index

    @staticmethod
    def __overlap_section(qmatoms, inp):
        pass

    def add(self, selection=None, is_link=False):
        qdf = Preparation.__clean_qdf(self.selector.select(selection))
        qdf.insert(2, 'is_link', [int(is_link)]*len(qdf))
        self.qm_atoms = self.qm_atoms.append(qdf)

    def delete(self, selection=None):
        qdf = Preparation.__clean_qdf(self.selector.select(selection))
        self.qm_atoms = self.qm_atoms.drop(qdf.index, errors='ignore')

    def clear(self, dummy=None):  # dummy makes main look nicer
        self.qm_atoms = pd.DataFrame()

    def view(self, dummy=None):  # dummy makes main look nicer
        print(self.qm_atoms)

    def prepare_cpmd(self, inp_tmp=None, mdp_inp=None, ndx_out=None, inp_out=None):
        """Args:
            inp_tmp: cpmd input file, used as template
            mdp_inp: gromacs input file, checked for errors
            ndx_out: gromacs index file, output
            inp_out: new mimic cpmd input file, output
        """
        # Check for obvious errors in selection
        if self.qm_atoms.empty:
            raise SelectionError('No atoms have been selected for the QM partition.')

        # Delete self.mpt for better garbage collection
        try:
            del self.mpt
        except AttributeError:
            pass

        # Retrieve number of steps and timestep from mdp_inp and do some checks
        try:
            mdp_inp = Mdp.from_file(mdp_inp)
            maxsteps, timestep, mdp_errors = mdp_inp.check()
        except:
            maxsteps, timestep, mdp_errors = 1000, 5.0, None
        if mdp_errors:
            logging.warning('Found md parameters which are inconsistent with MiMiC runs:')
            for error in mdp_errors:
                logging.warning(f'\t{error}')

        # Create an index group in GROMACS format (and write it to a file)
        qm_atoms_ndx_group = Preparation.__ndx_group(self.qm_atoms.index, 'QMatoms')
        if ndx_out:
            gbl.host.write(qm_atoms_ndx_group, ndx_out)
            logging.info(f'Wrote Gromacs index file to {ndx_out}')

        # index is also reset, for getOverlaps_Atoms()
        sorted_qm_atoms = self.qm_atoms.sort_values(by=['is_link', 'element']).reset_index()

        # Prepare CPMD input file
        try:
            cpmd = InputScript.from_file(inp_tmp)
        except:
            cpmd = InputScript()
        cpmd.mimic = Section()
        cpmd.mimic.paths = f'1\n{gbl.host.pwd()}'
        cpmd.mimic.box = '  '.join([str(s/BOHR_RADIUS) for s in self.structure.box])

#        cpmd = _qmhelper.getOverlaps_Atoms(sorted_qm_atoms, inp) # TODO: Move to scripts.cpmd

        cpmd.mimic.long_range__coupling = ''
        cpmd.mimic.fragment__sorting__atom_wise__update = 100
        cpmd.mimic.cutoff__distance = 20.0
        cpmd.mimic.multipole__order = 3

        total_charge = sum(self.qm_atoms['charge'])
        if not round(total_charge, 2).is_integer():
            logging.warning(f'Total charge of QM region is {total_charge}. Rounding to integer.')

        cpmd.system = Section()
        cpmd.system.charge = round(total_charge)

        # set cpmd section
        if not cpmd.checkSection('cpmd'):
            cpmd.cpmd = Section()
        cpmd.cpmd.mimic = ''
        cpmd.cpmd.parallel__constraints = ''

        if not cpmd.checkSection('dft'):
            cpmd.dft = Section()
        cpmd.dft.functional__blyp = '' # TODO: update to new XC_DRIVER code


        cpmd.cpmd.maxsteps = maxsteps
        cpmd.cpmd.timestep = round(timestep)

        if inp_out is None:
            logging.info('Created new CPMD input script for MiMiC run')
        else:
            gbl.host.write(str(cpmd), inp_out)
            logging.info(f'Wrote new CPMD input script to {inp_out}')

        return qm_atoms_ndx_group, cpmd
