import pandas as pd
from .selector import GroSelector
from ..io.mpt import Mpt
from ..io.gro import Gro
from ..scripts.mdp import Mdp
from ..scripts.cpmd import InputScript, Section
from .._global import _Global as gbl
from ..utils.errors import SelectionError
from ..utils.constants import BOHR_RADIUS, ATOMIC_TIME_UNIT


class Preparation:

    def __init__(self, topology, structure, selector=''):

        if isinstance(topology, Mpt):
            self.mpt = topology
        elif isinstance(topology, str):
            self.mpt = Mpt.from_file(topology)
        else:  # Raise exception
            pass

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
    def __check_mdp(mdp):
        nsteps = 1000
        dt = 5.0

        if mdp is None:
            return nsteps, dt, "\tNo mdp file given, using default values for number of steps and timestep"

        mdp_errors = []

        if not mdp.hasparam('integrator') or mdp.integrator != 'mimic':
            mdp_errors.append("\tWrong integrator for MiMiC run, set integrator = mimic")

        if not mdp.hasparam('qmmm_grps') or mdp.qmmm_grps != 'QMatoms':
            mdp_errors.append("\tIndex group for QM atoms does not correspond to ndx file, set QMMM-grps = QMatoms")

        if mdp.hasparam('constraints') and mdp.constraints != 'none':
            mdp_errors.append("\tMolecule should not be constrained, set constraints = none")

        if mdp.hasparam('tcoupl') and mdp.tcoupl != 'no':
            mdp_errors.append("\tTemperature coupling will not be active, set tcoupl = no")

        if mdp.hasparam('pcoupl') and mdp.pcoupl != 'no':
            mdp_errors.append("\tPressure coupling will not be active, set pcoupl = no")

        # TODO: Check for more errors in mdp file

        if mdp.hasparam('nsteps'):
            nsteps = int(mdp.nsteps)
        else:
            mdp_errors.append("\tNumber of steps is not given, using default value")

        if mdp.hasparam('dt'):
            dt = float(mdp.dt)/ATOMIC_TIME_UNIT
        else:
            mdp_errors.append("\tTimestep is not given, using default value")

        return nsteps, dt, "\n".join(mdp_errors)


    @staticmethod
    def __ndx_group(qm_ids, group_name):

        col_len = 15
        space_len = 6
        max_len = len(str(max(qm_ids)))
        spaces = space_len if max_len <= space_len else max_len

        index = f'[ {group_name} ]\n'
        for i, idx in enumerate(qm_ids):
            if i%col_len == 0:
                index += '\n'
            index += "{:{}}".format(idx, spaces)

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


    def clear(self):
        self.qm_atoms = pd.DataFrame()


    def prepare_cpmd(self, inp_tmp=None, mdp_inp=None, ndx_out=None, inp_out=None):
        """Args:
            inp_tmp: name of cpmd input file, used as template
        """
        # Check for obvious errors in selection
        if self.qm_atoms.empty:
            raise SelectionError("No atoms have been selected for the QM partition.")

        # Delete self.mpt for better garbage collection
        del self.mpt

        # If provided, retrieve number of steps and timestep from mdp_inp and do some checks
        if mdp_inp is not None:
            mdp_inp = Mdp.from_file(mdp_inp)
        maxsteps, timestep, mdp_errors = Preparation.__check_mdp(mdp_inp)

        if mdp_errors != "":
            gbl.logger.write('warning', '')
            gbl.logger.write('warning', f"Check your parameters in {mdp_inp}:")
            gbl.logger.write('warning', mdp_errors)
            gbl.logger.write('warning', '')

        # Create an index group in GROMACS format (and write it to a file)
        qm_atoms_ndx_group = Preparation.__ndx_group(self.qm_atoms.index, 'QMatoms')
        if ndx_out is not None:
            gbl.host.write(qm_atoms_ndx_group, ndx_out)
            gbl.logger.write('info', f"Wrote Gromacs index file to {ndx_out}")

        # index is also reset, for getOverlaps_Atoms()
        sorted_qm_atoms = self.qm_atoms.sort_values(by=['is_link', 'element']).reset_index()

        # Prepare CPMD input file
        if inp_tmp is not None:
            cpmd = InputScript.from_file(inp_tmp)
        else:
            cpmd = InputScript()

        cpmd.mimic = Section()
        cpmd.mimic.paths = f"1\n{gbl.host.pwd()}"
        cpmd.mimic.box = '  '.join([str(s/BOHR_RADIUS) for s in self.structure.box])

        inp = _qmhelper.getOverlaps_Atoms(sorted_qm, inp)

        inp.mimic.long_range__coupling = ''
        inp.mimic.fragment__sorting__atom_wise__update = 100
        inp.mimic.cutoff__distance = 20.0
        inp.mimic.multipole__order = 3

        q = sum(self.qm_atoms['charge'])
        if not round(q, 2).is_integer():
            gbl.logger.write('warning', (f'Total charge of QM region (={q}) not an integer up to 2 decimal places.'
                                         ' \nRounding to integer anyways..'))

        inp.system.charge = round(q) # system section already created in getQverlap_Atoms()

        # set cpmd section
        if not inp.checkSection('cpmd'):
            inp.cpmd = cpmd.Section()
        inp.cpmd.mimic = ''
        inp.cpmd.parallel__constraints = ''

        if not inp.checkSection('dft'):
            inp.dft = cpmd.Section()
        inp.dft.functional__blyp = '' # TO DO: update to new XC_DRIVER code

        # default values
        inp.cpmd.maxsteps = maxsteps
        inp.cpmd.timestep = round(timestep)

        if cpmd_file is None:
            gbl.logger.write('info', "Created CPMD input script..")
        else:
            gbl.host.write(str(inp), cpmd_file)
            gbl.logger.write('info', f"Wrote CPMD input script to {cpmd_file}..")

        return ndx, inp
