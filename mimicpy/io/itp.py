''' Module for itp files '''

import re
import pandas as pd
from os.path import dirname
from .parser import Parser
from ..utils.strs import clean
from ..utils.elements import elements


class Itp:
    ''' reads itp files '''

    columns = ['number', 'type', 'resid', 'resname', 'name', 'charge', 'element', 'mass']

    def __init__(self, file):
        self.file = file


    @staticmethod
    def section_is_in_string(section, string):
        if section == '*': # why?
            section = ''
        return bool('[' in string and ']' in string and section in string)


    @staticmethod
    def get_section(section, string, comments=';'):
        # clean string
        # find text b/w [ section ] and either [ or # or EOF
        # look for [ section ] / look for lines / look for optional spaces and [ or #
        string = clean(string, comments)
        section_regex = re.compile(fr"\[\s*{section}\s*\]\n((?:.+\n)+?)\s*(?:$|\[|#)", re.MULTILINE)
        section = section_regex.findall(string)
        return section


    @staticmethod
    def parse_block_till_section(itp, *sections, comments=[';', '#']):

        def get_section_headers(string):
            string = clean(string, comments)
            section_header_regex = re.compile(fr"\[\s*(.*?)\s*\]", re.MULTILINE)
            section_headers = section_header_regex.findall(string)
            return section_headers

        part_of_itp = ''
        read_sections = []

        for chunk in itp:
            part_of_itp += chunk
            sections_in_chunk = get_section_headers(chunk)
            read_sections += sections_in_chunk
            # all wanted section headers have been passed and last section has been read completely
            if set(sections).issubset(read_sections) and read_sections[-1] != sections[-1]:
                break

        return part_of_itp


    def load_molecules_and_atoms(self, buffer=1000):
        itp_file = Parser(self.file, buffer)
        itp_text = ''

        while not itp_file.is_closed:
            block = self.parse_block_till_section(itp_file, 'moleculetype', 'atoms')
            if any([self.section_is_in_string(section, block) for section in ['moleculetype', 'atoms']]):
                itp_text += block

        return itp_text


    def get_included_topology_files(self, string, comments=';'):
        string = clean(string, comments)
        include_file_regex = re.compile(r"#include\s+[\"\'](.+)\s*[\"\']", re.MULTILINE)
        included_itps = include_file_regex.findall(string)
        included_itps = [(dirname(self.file) + '/' + itp) for itp in included_itps]
        return included_itps


    def get_atomtypes(self, buffer=1000):
        itp_file = Parser(self.file, buffer)
        itp_text = self.parse_block_till_section(itp_file, 'atomtypes')
        atomtypes_section = self.get_section('atomtypes', itp_text)

        if atomtypes_section == []: # What if atomtypes are in several itps?
            included_itps = self.get_included_topology_files(itp_text)

            for included_itp in included_itps:
                itp = Itp(included_itp)
                atom_types = itp.get_atomtypes(buffer=buffer)
                if atom_types != {}:
                    return atom_types

            return {}
        else:
            atomtypes_section = atomtypes_section[0]

        atom_types = {}

        for line in atomtypes_section.splitlines():
            line = line.split()
            atom_type = line[0]
            atom_number = line[1]
            if not atom_number.isnumeric():
                # raise Exception
                continue
            atom_types[atom_type] = elements[int(atom_number)]

        return atom_types


    def read(self, requested_molecules, atom_types, buffer=1000):

        def guess_element_from(mass, name, atom_type): # Make a huge fuss about guessing elements
                mass_int = int(mass)
                
                if mass_int <= 0:
                    if name in elements:
                        elem = name
                    elif  atom_type in elements:
                        atom_type = name
                if mass_int<=1: elem = 'H' # for H
                elif mass_int<36: elem = elements[mass_int//2] # He to Cl
                else: elem = name.title() # from Ar onwards, assume name same as symbol, case insensitive
                return elem

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
                else: # Give warning
                    pass

                number = int(number)
                resid = int(resid)
                charge = float(charge)
                mass = float(mass)

                if atom_type in atom_types:
                    element = atom_types[atom_type]
                else:
                    element = guess_element_from(mass, name, atom_type )

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

        itp_text = self.load_molecules_and_atoms(buffer)
        molecule_section = self.get_section('moleculetype', itp_text)
        atom_section = self.get_section('atoms', itp_text, comments=[';', '#'])

        if molecule_section == [] and atom_section == []:
            return None

        molecules = []
        atom_infos = []

        for molecule, atoms in zip(molecule_section, atom_section):

            mol = molecule.split()[0]
            if mol not in requested_molecules:
                continue
            molecules.append(mol)
            atom_infos.append(read_atoms(atoms))        

        return dict(zip(molecules, atom_infos))
