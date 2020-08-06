""" Module for itp files """

from os.path import dirname
import re
import pandas as pd
from .parser import Parser
from ..utils.strs import clean
from ..utils.elements import elements


class Itp:
    """ reads itp files """

    columns = ['number', 'type', 'resid', 'resname', 'name', 'charge', 'element', 'mass']

    def __init__(self, file, requested_molecules=None, atom_types=None, buffer=1000, mode='r'):
        self.file = file
        self.requested_molecules = requested_molecules
        self.atom_types_dict = atom_types
        self.buffer = 1000
        self.mode = mode
        
        self._topol = None
        
        self._topology_files = None
        self._molecules = None
        
        if mode == 'r':
            self.__read()
        elif mode == 'w':
            pass
        elif mode == 't': # as topol.top
            self.__read_as_topol()
        else:  # Raise exception
            pass
    
    ### Static Helper Methods
    #
    @staticmethod
    def __get_molecules(topology):
        molecules = []

        for line in topology.splitlines()[::-1]:
            if Itp.__section_is_in_string('molecules', line):
                break
            if line.strip() == '' or line.startswith(';'):
                continue
            molecule_name, number_of_molecules = line.split()
            molecules += [(molecule_name, int(number_of_molecules))]

        return molecules[::-1]

    @staticmethod
    def __section_is_in_string(section, string):
        return bool('[' in string and ']' in string and section in string)


    @staticmethod
    def __get_section(section, string, comments=';'):
        # Clean string
        # Find text b/w [ section ] and either [ or # or EOF
        # Look for [ section ] / look for lines / look for optional spaces and [ or #
        string = clean(string, comments)
        section_regex = re.compile(fr"\[\s*{section}\s*\]\n((?:.+\n)+?)\s*(?:$|\[|#)", re.MULTILINE)
        section = section_regex.findall(string)
        return section


    @staticmethod
    def __parse_block_till_section(itp, *sections, comments=[';', '#']):

        def __get_section_headers(string):
            string = clean(string, comments)
            section_header_regex = re.compile(fr"\[\s*(.*?)\s*\]", re.MULTILINE)
            section_headers = section_header_regex.findall(string)
            return section_headers

        part_of_itp = ''
        read_sections = []

        for chunk in itp:
            part_of_itp += chunk
            sections_in_chunk = __get_section_headers(chunk)
            read_sections += sections_in_chunk
            # All wanted section headers have been passed and last section has been read completely
            if set(sections).issubset(read_sections) and read_sections[-1] != sections[-1]:
                break

        return part_of_itp
    
    ##Non static Helper Methods -- These use self.file and/or self.buffer
    #
    def __load_molecules_and_atoms(self):
        itp_file = Parser(self.file, self.buffer)
        itp_text = ''

        while not itp_file.is_closed:
            block = self.__parse_block_till_section(itp_file, 'moleculetype', 'atoms')
            if any([self.__section_is_in_string(section, block) for section in ['moleculetype', 'atoms']]):
                itp_text += block

        return itp_text

    def __get_included_topology_files(self, string, comments=';'):
        string = clean(string, comments)
        include_file_regex = re.compile(r"#include\s+[\"\'](.+)\s*[\"\']", re.MULTILINE)
        included_itps = include_file_regex.findall(string)
        included_itps = [(dirname(self.file) + '/' + itp) for itp in included_itps]
        return included_itps
    
    ### Read functions
    #
    def __read_atomtypes(self):
        itp_file = Parser(self.file, self.buffer)
        itp_text = self.__parse_block_till_section(itp_file, 'atomtypes')
        atomtypes_section = self.__get_section('atomtypes', itp_text)

        if atomtypes_section == []: # What if atomtypes are in several itps?
            included_itps = self.__get_included_topology_files(itp_text)

            for included_itp in included_itps:
                itp = Itp(included_itp)
                atom_types = itp.__read_atomtypes()
                if atom_types != {}:
                    return atom_types

            return {}

        atomtypes_section = atomtypes_section[0]

        atom_types = {}

        for line in atomtypes_section.splitlines():
            line = line.split()
            atom_type = line[0]
            atom_number = line[1]
            if not atom_number.isnumeric():
                # Raise exception
                continue
            atom_types[atom_type] = elements[int(atom_number)]

        return atom_types
    
    def __read(self):
        
        if self.atom_types_dict is None or self.requested_molecules is None:
            pass # raise MiMiCPyError

        def guess_element_from(mass, name, atom_type):  # Make a huge fuss about guessing elements
            mass_int = int(mass)

            if mass_int <= 0:  # Cannot guess from mass
                if name in elements.values():  # Guess from atom name
                    element = name
                elif atom_type in elements.values():  # Guess from atom type
                    element = atom_type

            elif mass_int <= 1:  # Guess H from mass
                element = 'H'

            elif mass_int < 36:  # Guess He to Cl from mass
                element = elements[mass_int//2]

            else:  # Assume name is element symbol from Ar onwards
                element = name.title()  # Case insensitive

            return element

        def read_atoms(atom_section):
            cols = self.columns
            atom_info = {k:[] for k in cols}

            for _, line in enumerate(atom_section.splitlines()):
                line = line.split()
                if len(line) == 8:
                    number, atom_type, resid, resname, name, _, charge, mass = line[:8]
                elif len(line) == 7:
                    number, atom_type, resid, resname, name, _, charge = line[:7]
                    mass = 0
                else:  # Give warning
                    pass

                number = int(number)
                resid = int(resid)
                charge = float(charge)
                mass = float(mass)

                if atom_type in self.atom_types_dict:
                    element = self.atom_types_dict[atom_type]
                else:
                    element = guess_element_from(mass, name, atom_type)

                atom_info[cols[0]].append(number)
                atom_info[cols[1]].append(atom_type)
                atom_info[cols[2]].append(resid)
                atom_info[cols[3]].append(resname)
                atom_info[cols[4]].append(name)
                atom_info[cols[5]].append(charge)
                atom_info[cols[6]].append(element)
                atom_info[cols[7]].append(mass)

            atoms = pd.DataFrame(atom_info).set_index(cols[0])
            return atoms

        itp_text = self.__load_molecules_and_atoms()
        molecule_section = self.__get_section('moleculetype', itp_text)
        atom_section = self.__get_section('atoms', itp_text, comments=[';', '#'])

        if molecule_section == [] and atom_section == []:
            return None

        molecules = []
        atom_infos = []

        for molecule, atoms in zip(molecule_section, atom_section):

            mol = molecule.split()[0]
            if mol not in self.requested_molecules:
                continue
            molecules.append(mol)
            atom_infos.append(read_atoms(atoms))

        self._topol = dict(zip(molecules, atom_infos))
    
    def __read_as_topol(self):
        top_parser = Parser(self.file)
        topology = ''.join(top_parser)

        self._molecules = self.__get_molecules(topology)

        self._topology_files = self.__get_included_topology_files(topology)
        self._topology_files.append(self.file)

    ###Property getters
    #
    @property
    def topol(self):
        if self.mode == 'r':
            return self._topol
        
        self.mode = 'r'
        self.__read()
        return self._topol
    
    @property
    def molecules(self):
        if self.mode == 't':
            return self._molecules
        
        self.mode = 't'
        self.__read_as_topol()
        return self._molecules
    
    @property
    def topology_files(self):
        if self.mode == 't':
            return self._topology_files
        
        self.mode = 't'
        self.__read_as_topol()
        return self._topology_files
    
    @property
    def atom_types(self):
        return self.__read_atomtypes()
    #
    ###
