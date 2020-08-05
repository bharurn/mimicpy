""" Module for top files """

from .parser import Parser
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
        """ reads molecule and atom information """

        def get_molecules(topology):
            molecules = []

            for line in topology.splitlines()[::-1]:
                if Itp.section_is_in_string('molecules', line):
                    break
                if line.strip() == '' or line.startswith(';'):
                    continue
                molecule_name, number_of_molecules = line.split()
                molecules += [(molecule_name, int(number_of_molecules))]

            return molecules[::-1]

        top_parser = Parser(self.file)
        topology = ''.join(top_parser)

        molecules = get_molecules(topology)
        molecule_types = [m[0] for m in molecules]

        top = Itp(self.file)
        topology_files = top.get_included_topology_files(topology)
        topology_files.append(self.file)

        itp = Itp(topology_files[0])
        atom_types = itp.get_atomtypes()
        if self.nonstandard_atom_types is not None:
            atom_types.update(self.nonstandard_atom_types)

        atoms = {}

        for itp in topology_files:
            itp = Itp(itp)
            atom_info = itp.read(molecule_types, atom_types, self.buffer)
            if atom_info:
                atoms.update(atom_info)

        topol_dict = TopolDict.from_dict(atoms)

        self._molecules = molecules
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
