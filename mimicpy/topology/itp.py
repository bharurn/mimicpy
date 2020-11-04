"""Module for itp files"""

from os.path import dirname, isfile, join
import re
import logging
import pandas as pd
from ..utils.file_handler import Parser
from ..utils.strings import clean
from ..utils.elements import ELEMENTS
from ..utils.errors import MiMiCPyError, ParserError

class Itp:
    """reads itp files"""

    columns = ['number', 'type', 'resid', 'resname', 'name', 'charge', 'element', 'mass']

    def __init__(self, file, requested_molecules=None, atom_types=None, buffer=1000, mode='r', guess_elements=True, gmxdata=''):
        self.file = file
        self.requested_molecules = requested_molecules
        self.atom_types_dict = atom_types
        self.buffer = buffer
        self.mode = mode
        self.guess_elements = guess_elements
        self.gmxdata = gmxdata
        self._topol = None
        self._topology_files = None
        self._molecules = None
        self._molecule_types = None
        self.guessed_elems_history = {}

        if mode == 'r':
            self.__read()
        elif mode == 't':
            self.__read_as_topol()
        elif mode == 'w':
            pass
        else:
            raise MiMiCPyError(f'{mode} is not a mode. Only r, t, or w can be used.')

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
    def molecule_types(self):
        if self.mode == 't':
            return self._molecule_types
        self.mode = 't'
        self.__read_as_topol()
        return self._molecule_types

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
    def __get_section(section, string):
        # Clean string - or not?
        # Find text b/w [ section ] and either [ or # or EOF
        # Look for [ section ] / look for lines / look for optional spaces and [ or #
        # string = clean(string, comments)
        section_regex = re.compile(fr"\[\s*{section}\s*\]\n((?:.+\n)+?)\s*(?:$|\[|#)", re.MULTILINE)
        section_list = section_regex.findall(string)
        return section_list

    @staticmethod
    def __parse_block_till_section(itp, *sections):

        def get_section_headers(string):
            section_header_regex = re.compile(fr"\[\s*(.*?)\s*\]", re.MULTILINE)
            section_headers = section_header_regex.findall(string)
            return section_headers

        part_of_itp = ''
        read_sections = []
        for chunk in itp:
            part_of_itp += chunk
            sections_in_chunk = get_section_headers(chunk)
            read_sections += sections_in_chunk
            # All wanted section headers have been passed and last section has been read completely
            if set(sections).issubset(read_sections) and read_sections[-1] != sections[-1]:
                break
        return part_of_itp

    def __load_molecules_and_atoms(self):
        itp_file = Parser(self.file, self.buffer)
        itp_text = ''
        while not itp_file.is_closed:
            block = Itp.__parse_block_till_section(itp_file, 'moleculetype', 'atoms')
            if any([Itp.__section_is_in_string(section, block) for section in ['moleculetype', 'atoms']]):
                itp_text += block
        return itp_text

    def __get_included_topology_files(self, string, comments=';'):
        string = clean(string, comments)
        include_file_regex = re.compile(r"#include\s+[\"\'](.+)\s*[\"\']", re.MULTILINE)
        included_itps = include_file_regex.findall(string)

        if self.gmxdata is None:
            included_itps = [join(dirname(self.file), itp) for itp in included_itps]
        else:
            included_itps = [join(dirname(self.file), itp) if isfile(join(dirname(self.file), itp)) \
                         else join(self.gmxdata, itp) for itp in included_itps]
        return included_itps
    
    def __get_all_atomtypes_sections(self):
        itp_file = Parser(self.file, self.buffer)
        itp_text = Itp.__parse_block_till_section(itp_file, 'atomtypes')
        clean_itp_text = clean(itp_text, comments=';')
        atomtypes_section = "\n".join(Itp.__get_section('atomtypes', clean_itp_text))
        if atomtypes_section == "":
            included_itps = self.__get_included_topology_files(clean_itp_text)
            for included_itp in included_itps:
                try:
                    itp = Itp(included_itp)
                    atom_types = itp.__get_all_atomtypes_sections()
                    if atom_types is not None:
                        atomtypes_section += atom_types
                except OSError:
                    logging.warning('Could not find %s. Skipping.', included_itp)
            return atomtypes_section
        else:
           return atomtypes_section


    def __read_atomtypes(self):
        atomtypes_section = self.__get_all_atomtypes_sections()
        cols = ['type', 'X', 'mass', 'charge', 'ptype', 'sigma', 'epsilon']
        float_cols = [cols[i] for i in range(7) if i in [2,3,5,6]]
        atm_types_section_dct = {k:[] for k in cols}
        for line in atomtypes_section.splitlines():
            line_split = line.split()
            if len(line_split) != 7:
                raise Exception
            else:
                for i, col in enumerate(cols):
                    atm_types_section_dct[col].append(float(line_split[i]) if col in float_cols else line_split[i])
        
        self.atom_types_df = pd.DataFrame(atm_types_section_dct)
        df = self.atom_types_df.copy()
        
        convert_elem = lambda x: ELEMENTS[int(x)] if x.isnumeric() and int(x) > 0 else None
        elem_col = df.columns[1]
        
        df[elem_col] = df[elem_col].apply(convert_elem)
        df = df[df[elem_col].notnull()]
        df = df.set_index(cols[0])
        return df[elem_col].to_dict()

    def __read(self):

        def guess_element_from(mass, name, atom_type):
            element = 'H'
            mass_int = int(round(mass))
            if mass_int <= 0:  # Cannot guess from mass
                if name in ELEMENTS.values():  # Guess from atom name
                    element = name
                elif atom_type in ELEMENTS.values():  # Guess from atom type
                    element = atom_type
                else:
                    element = name.title()[0]
            elif mass_int <= 1:  # Guess H from mass
                element = 'H'
            elif mass_int < 36:  # Guess He to Cl from mass
                element = ELEMENTS[mass_int//2]
            else:  # Assume name is element symbol from Ar onwards
                element = name.title()  # Case insensitive
            self.guessed_elems_history[atom_type] = element
            return element

        def read_atoms(atom_section):
            cols = self.columns
            atom_info = {k:[] for k in cols}
            number_of_bad_lines = 0
            for i, line in enumerate(atom_section.splitlines()):
                line = line.split()
                if len(line) == 8:
                    number, atom_type, resid, resname, name, _, charge, mass = line[:8]
                elif len(line) == 7:
                    number, atom_type, resid, resname, name, _, charge = line[:7]
                    mass = 0
                else:
                    if number_of_bad_lines > 5:
                        raise ParserError(self.file, 'topology',
                                          details='Too many bad lines in [ atoms ] section.')
                    logging.warning('Atom %s in [ atoms ] section is not formatted properly. Skipping...', i)
                    number_of_bad_lines += 1
                    continue
                number = int(number)
                resid = int(resid)
                charge = float(charge)
                mass = float(mass)
                if self.atom_types_dict is not None and atom_type in self.atom_types_dict:
                    element = self.atom_types_dict[atom_type]
                elif self.guess_elements:
                    element = guess_element_from(mass, name, atom_type)
                else:
                    raise ParserError(self.file, 'topology', details=('Cannot determine atomic symbol'
                                     f' for atom with name {name} and type {atom_type} in residue {resname}'))
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
        clean_itp_text = clean(itp_text,  comments=[';', '#'])
        molecule_section = Itp.__get_section('moleculetype', clean_itp_text)
        atom_section = Itp.__get_section('atoms', clean_itp_text)
        if molecule_section == [] and atom_section == []:
            return None
        molecules = []
        atom_infos = []

        for molecule, atoms in zip(molecule_section, atom_section):
            mol = molecule.split()[0]
            if self.requested_molecules is not None and mol not in self.requested_molecules:
                continue
            molecules.append(mol)
            atom_infos.append(read_atoms(atoms))
        self._topol = dict(zip(molecules, atom_infos))

    def __read_as_topol(self):
        top_parser = Parser(self.file)
        topology = ''.join(top_parser)
        self._molecules = Itp.__get_molecules(topology)
        self._molecule_types = [m[0] for m in self._molecules]
        self._topology_files = self.__get_included_topology_files(topology)
        self._topology_files.append(self.file)
        self.requested_molecules = self._molecule_types
