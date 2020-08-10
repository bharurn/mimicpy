""" Module for top files """

from .itp import Itp
from .topol_dict import TopolDict


class Top:
    """ reads top files """

    def __init__(self, file, mode='r', buffer=1000, nonstandard_atom_types=None):
        self.file = file
        self.mode = mode
        self.buffer = buffer
        self.nonstandard_atom_types = nonstandard_atom_types
        self._molecules = None
        self._topol_dict = None

        if mode == 'r':
            self.__read()
        elif mode == 'w':
            pass
        else:  # Raise exception
            pass


    def __read(self):
        """ Read molecule and atom information """

        top = Itp(self.file, mode='t')
        atom_types = top.atom_types

        if self.nonstandard_atom_types is not None:
            atom_types.update(self.nonstandard_atom_types)

        atoms = {}
        molecule_types = top.molecule_types

        for itp in top.topology_files:
            itp = Itp(itp, molecule_types, atom_types, self.buffer)
            if itp.topol is not None:
                atoms.update(itp.topol)

        topol_dict = TopolDict.from_dict(atoms)

        self._molecules = top.molecules
        self._topol_dict = topol_dict


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
