from abc import ABC, abstractmethod
from ..utils.errors import MiMiCPyError, ParserError
from ..utils.file_handler import Parser, write as write_string

class BaseCoordsClass(ABC):
    def __init__(self, file_name, buffer=1000):
        self.file_name = file_name
        self.buffer = buffer
    
    def read(self):
        self.file = Parser(self.file_name)
        coords, box = self._read()
        self.file.close()
        return coords, box
        
    @abstractmethod
    def _read(self):
        pass

    def write(self, sele, coords=None, box=None):
        if coords:
            sele = sele.join(coords)
         
        write_string(self._write(sele.reset_index(), box), self.file_name, 'w')

    @abstractmethod
    def _write(self, mpt_coords, box):
        pass

# adapter class
class CoordsIO:
    def __init__(self, file_name, mode='r', buffer=1000, ext=None):
        if isinstance(file_name, BaseCoordsClass):
            self.__coords_obj = file_name
        else:
            if ext is None:
                ext = file_name.split('.')[-1]
            
            if ext == 'gro':
                from .gro import Gro
                self.__coords_obj = Gro(file_name, buffer)
            elif ext == 'pdb':
                from .pdb import Pdb
                self.__coords_obj = Pdb(file_name, buffer)
            else:
                raise ParserError('Unknown coordinate format')
                
        self.mode = mode
        self._coords = None
        self._box = None

        if mode == 'r':
            self.__read()
        elif mode == 'w':
            pass
        else:
            raise MiMiCPyError(f'{mode} is not a mode. Only r or w can be used.')

    @property
    def coords(self):
        if self.mode != 'r': self.__read()
        return self._coords

    @property
    def box(self):
        if self.mode != 'r': self.__read()
        return self._box

    def __read(self):
        self.mode = 'r'
        self._coords, self._box = self.__coords_obj.read()

    def write(self, sele, coords=None, box=None):
        if self.mode != 'w':
            self.mode = 'w'
        self.__coords_obj.write(sele, coords, box)