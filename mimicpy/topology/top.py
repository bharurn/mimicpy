"""Module for top files"""

import logging
from os import environ
from os.path import basename, join
from .itp import Itp
from .topol_dict import TopolDict
from ..utils.errors import MiMiCPyError


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

        top = Itp(self.file, mode='t', gmxdata=self.gmxdata)
        atom_types = top.atom_types
        molecule_types = top.molecule_types

        if self.nonstandard_atomtypes is not None:
            # TODO: Support non-standard atomtypes input as file (list, itp, ...)
            atom_types.update(self.nonstandard_atomtypes)

        atoms = {}

        for itp_file in top.topology_files:
            itp_file_name = basename(itp_file) # print only file name, and not full path
            try:
                itp = Itp(itp_file, molecule_types, atom_types, self.buffer, 'r', self.guess_elements, self.gmxdata)
                if itp.topol is not None:
                    atoms.update(itp.topol)
                    logging.info('Read atoms from %s.', itp_file_name)
                else:
                    logging.info('No atoms found in %s.', itp_file_name)
            except OSError:
                logging.warning('Could not find %s in local or Gromacs data directory. Skipping...', itp_file_name)
        topol_dict = TopolDict.from_dict(atoms)

        self._molecules = top.molecules
        self._topol_dict = topol_dict