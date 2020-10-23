from importlib import import_module
from abc import ABC, abstractmethod
from ..utils.errors import MiMiCPyError
from ..topology.mpt import Mpt
from ..utils.file_handler import Parser, write as write_string

class BaseCoordsClass(ABC):
    def __init__(self, file_name, buffer):
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

    def write(self, coords, mpt_file, ids=None):
        if ids is None:
            ids = coords['ids']
        else:
            coords['id'] = ids

        mpt = Mpt.fromFile(mpt_file)
        mpt_coords = mpt[ids].merge(coords, left_on='id', right_on='id').set_index(['id'])

        write_string(self._write(mpt_coords), self.file_name, 'w')

    @abstractmethod
    def _write(self, mpt_coords):
        pass

# adapter class
class CoordsIO(ABC):
    def __init__(self, file_name, mode='r', buffer=1000, ext=None):
        if ext is None:
            ext = file_name.split('.')[-1]

        self.coord_obj = import_module('.'+ext+'.'+ext.title())(file_name, mode, buffer)
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
        self._coords, self._box = self.coords_obj.read()

    def write(self, coords, mpt_file, ids=None):
        if self.mode != 'w':
            self.mode = 'w'
        self.coords_obj.write(coords, mpt_file, ids)