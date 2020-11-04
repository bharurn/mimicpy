"""Module for top files"""

import logging
from os import environ
from os.path import basename, join
from .itp import Itp
from .topol_dict import TopolDict
from ..utils.errors import MiMiCPyError
from ..utils.strings import print_dict
from ..utils.atomic_numbers import atomic_numbers
from ..utils.file_handler import write


class Top:
    """reads top files"""

    def __init__(self, file, mode='r', buffer=1000, nonstandard_atomtypes=None, guess_elements=True, gmxdata=None):
        self.file = file
        self.mode = mode
        self.buffer = buffer
        self.nonstandard_atomtypes = nonstandard_atomtypes
        self.guess_elements = guess_elements
        
        if gmxdata is None:
            if 'GMXDATA' in environ:
                gmxdata = join(environ['GMXDATA'], 'top')
            elif 'GMXLIB' in environ:
                gmxdata = join(environ['GMXLIB'], 'top')
        
        self.gmxdata = gmxdata
        
        if self.gmxdata: 
            logging.info('Using {} as path to Gromacs installation.'.format(self.gmxdata))
        else:    
            self.gmxdata = ''
            logging.warning('Cannot find path to Gromacs installation.')
            
        
        self._molecules = None
        self._topol_dict = None
        
        if mode == 'r':
            self.__read()
        elif mode == 'w':
            self.__read(True)
        else:
            raise MiMiCPyError(f'{mode} is not a mode. Only r or w can be used.')

    @property
    def molecules(self):
        if self.mode in ['r', 'w']:
            return self._molecules
        self.mode = 'r'
        self.__read()
        return self._molecules

    @property
    def topol_dict(self):
        if self.mode in ['r', 'w']:
            return self._topol_dict
        self.mode = 'r'
        self.__read()
        return self._topol_dict

    def __read(self, get_atomtypes=False):
        """Read molecule and atom information"""

        top = Itp(self.file, mode='t', gmxdata=self.gmxdata)
        atom_types = top.atom_types
        if get_atomtypes:
            self.atomtypes = top.atom_types_df
        else:
            self.atomtypes = None
        molecule_types = top.molecule_types

        if self.nonstandard_atomtypes is not None:
            atom_types.update(self.nonstandard_atomtypes)

        atoms = {}
        guessed_elems_history = {}
        
        for itp_file in top.topology_files:
            itp_file_name = basename(itp_file) # print only file name, and not full path
            try:
                itp = Itp(itp_file, molecule_types, atom_types, self.buffer, 'r', self.guess_elements, self.gmxdata)
                if itp.topol is not None:
                    atoms.update(itp.topol)
                    guessed_elems_history.update(itp.guessed_elems_history)
                    logging.info('Read atoms from %s.', itp_file_name)
                else:
                    logging.info('No atoms found in %s.', itp_file_name)
            except OSError:
                logging.warning('Could not find %s in local or Gromacs data directory. Skipping...', itp_file_name)
        topol_dict = TopolDict.from_dict(atoms)

        self._molecules = top.molecules
        self._topol_dict = topol_dict
        
        if guessed_elems_history:
            logging.warning('\nSome atom types had no atom numbers infomation.\nThey were guessed as follows:\n')
            print_dict(guessed_elems_history, "Atom Type", "Element", logging.warning)
    
    def write_atomtypes(self, file):
        if self.mode != 'w':
            return self.__read(True)
        
        elements = {}
        for k, df in self.topol_dict.todict().items():
            elements.update(dict(zip(df['type'], df['element'])))
        elements = {k:atomic_numbers[v] for k,v in elements.items()}
        
        itp_str = "[ atomtype ]\n"
        itp_str += ";  {:^11}{:^6}{:^10}{:^10}{:^6}     {}     {}\n".format('name','at.num','mass','charge','ptype',
                                                                           'sigma','epsilon')
        for i, row in self.atomtypes.iterrows():
            lst = [row[c] for c in self.atomtypes.columns]
            if lst[0] in elements:
                lst[1] = elements[lst[0]]
            else:
                lst[1] = int(lst[1])
            itp_str += "{:>11}{:6d}{:11.4f}{:11.4f}{:>6}     {:e}     {:e}\n".format(*lst)
            
        write(itp_str, file)