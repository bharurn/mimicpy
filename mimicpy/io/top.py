''' Module for top files '''

from .itp import Itp
from .topol_dict import TopolDict

class Top:
    ''' reads top files '''

    def __init__(self, file, mode='r'):
        self.file = file
        if mode == 'r':
            self.coords, self.box = self.__read()
        elif mode == 'w':
            pass
        else: # Raise Exception
            pass


    def read(self, buffer=1000, nonstandard_atom_types=None):
        ''' reads molecule and atom information '''

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

        with open(self.file, 'r') as f:
            topology = f.read()

        molecules = get_molecules(topology)
        molecule_types = [m[0] for m in molecules]

        top = Itp(self.file)
        topology_files = top.get_included_topology_files(topology)
        topology_files.append(self.file)

        itp = Itp(topology_files[0])
        atom_types = itp.get_atomtypes()
        if nonstandard_atom_types is not None:
            atom_types.update(nonstandard_atom_types)

        atoms = {}

        for itp in topology_files:
            itp = Itp(itp)
            atom_info = itp.read(molecule_types, atom_types, buffer)
            if atom_info:
                atoms.update(atom_info)

        topol_dict = TopolDict.from_dict(atoms)

        return molecules, topol_dict
