"""Module for top files"""

import logging
from .itp import Itp
from .topol_dict import TopolDict
from ..utils.errors import MiMiCPyError


class Top:
    """reads top files"""

    def __init__(self, file, mode='r', buffer=1000, nonstandard_atomtypes=None):
        self.file = file
        self.mode = mode
        self.buffer = buffer
        self.nonstandard_atomtypes = nonstandard_atomtypes
        self._molecules = None
        self._topol_dict = None

        if mode == 'r':
            self.__read()
        elif mode == 'w':
            pass
        else:
            raise MiMiCPyError(f'{mode} is not a mode. Only r or w can be used.')

    @property
    def molecules(self):
        if self.mode == 'r':
            return self._molecules
        self.mode = 'r'
        self.__read()
        return self._molecules

    @property
    def topol_dict(self):
        if self.mode == 'r':
            return self._topol_dict
        self.mode = 'r'
        self.__read()
        return self._topol_dict

    def __read(self):
        """Read molecule and atom information"""

        top = Itp(self.file, mode='t')
        atom_types = top.atom_types
        molecule_types = top.molecule_types

        if self.nonstandard_atomtypes is not None:
            # TODO: Support non-standard atomtypes input as file (list, itp, ...)
            atom_types.update(self.nonstandard_atomtypes)

        atoms = {}
        for itp in top.topology_files:
            try:
                itp = Itp(itp, molecule_types, atom_types, self.buffer)
                if itp.topol is not None:
                    atoms.update(itp.topol)
            except OSError:
                logging.warning('Could not find %s. Skipping.', itp)
        topol_dict = TopolDict.from_dict(atoms)

        self._molecules = top.molecules
        self._topol_dict = topol_dict
