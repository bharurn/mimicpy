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

    def get_mimic_input(self, inp_tmp=None, ndx_out=None, inp_out=None):
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
            raise SelectionError('No atoms have been selected for the QM partition')
            
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
        elif isinstance(inp_tmp, str):
            cpmd = CpmdScript.from_file(inp_tmp)
        else:
            cpmd = inp_tmp
        
        cpmd.atoms.clear_parameters() # clear atoms from inp_temp

        # Get overlaps and atoms
        overlaps = '{}'.format(len(sorted_qm_atoms))
        for i, atom in sorted_qm_atoms.iterrows():
            gromacs_id = atom['id']
            cpmd_id = i + 1
            overlaps += '\n2 {} 1 {}'.format(gromacs_id, cpmd_id)
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
            logging.warning('Total charge of QM region is %s, Rounding to integer', total_charge)
        cpmd.system.charge = round(total_charge)

        if not cpmd.cpmd.has_parameter('maxstep'):
            cpmd.cpmd.maxstep = 1000
        
        if not cpmd.cpmd.has_parameter('timestep'):
            cpmd.cpmd.timestep = 5.0

        if inp_out is None:
            logging.info('Created new CPMD input script for MiMiC run')
        else:
            write(str(cpmd), inp_out, 'w')
            logging.info('Wrote new CPMD input script to %s', inp_out)

        return qm_ndx_group, cpmd
    
    @staticmethod
    def get_gmx_input(inp=None, qmatoms=None, out=None):
        
        if qmatoms is None:
            qmatoms = 'qmatoms'
            
        if inp is None:
            mdp = Mdp()
        elif isinstance(inp, str):
            mdp = Mdp.from_file(inp)
        else:
            mdp = inp
            
        errors = False
            
        # TODO: Check for more errors in mdp file
        if (not mdp.has_parameter('integrator') or mdp.integrator != 'mimic') and inp != None:
            logging.warning('Wrong integrator for MiMiC run, setting integrator = mimic')
            errors = True
        
        if (not mdp.has_parameter('qmmm_grps') or mdp.qmmm_grps != qmatoms) and inp != None:
            logging.warning('Index group for QM atoms is not qmatoms, setting QMMM-grps to the appropriate group')
            errors = True
        
        if mdp.has_parameter('constraints') and mdp.constraints != 'none':
            logging.warning('Molecules should not be constrained by GROMACS, setting constraints = none')
            errors = True
            
        if mdp.has_parameter('tcoupl') and mdp.tcoupl != 'no':
            logging.warning('Temperature coupling will not be active, setting tcoupl = no')
            errors = True
            
        if mdp.has_parameter('pcoupl') and mdp.pcoupl != 'no':
            logging.warning('Pressure coupling will not be active, setting pcoupl = no')
            errors = True
            
        if mdp.has_parameter('dt'):
            dt = float(mdp.dt)
            if dt > 0.0001:
                logging.warning('Timestep may be too high for Gromacs (dt > 0.001)')
                errors = True
                
        mdp.integrator = 'mimic'
        mdp.dt = 0.0001
        mdp.constraints = 'none'
        mdp.tcoupl = 'no'
        mdp.pcoupl = 'no'
        mdp.qmmm_grps = qmatoms
        
        if not errors and inp != None:
            if isinstance(inp, str):
                fname = inp
            else:
                fname = 'MDP script'
            logging.info('No errors found in {}'.format(fname))
        elif out is None:
            logging.info('Created new MDP script for MiMiC run')
        else:
            write(str(mdp), out, 'w')
            logging.info('Wrote fixed MDP script to %s', out)