from importlib import import_module
from abc import ABC, abstractmethod
from ..utils.errors import MiMiCPyError
from ..topology.mpt import Mpt
from ..utils.file_handler import Parser, write

class Coords(ABC):
    def __init__(self, file, mode='r', buffer=1000):
        self.file_name = file
        self.mode = mode
        self.buffer = buffer
        self._coords = None
        self._box = None

        if mode == 'r':
            self.file = Parser(self.file_name)
            self._read()
        elif mode == 'w':
            pass
        else:
            raise MiMiCPyError(f'{mode} is not a mode. Only r or w can be used.')

    @property
    def coords(self):
        if self.mode == 'r':
            return self._coords
        self.mode = 'r'
        self._read()
        return self._coords

    @property
    def box(self):
        if self.mode == 'r':
            return self._box
        self.mode = 'r'
        self._read()
        return self._box
    
    @abstractmethod
    def _read(self):
        pass
    
    def write(self, file_name, coords, mpt_file, ids=None):
        if ids is None:
            ids = coords['ids']
        else:
            coords['id'] = ids
        
        mpt = Mpt.fromFile(mpt_file)
        mpt_coords = mpt[ids].merge(coords, left_on='id', right_on='id').set_index(['id'])
        
        write(self._write(mpt_coords), file_name, 'w')
    
    @abstractmethod
    def _write(self):
        pass

def read_coords(file, ext=None, mode='r', buffer=1000):
    if ext is None:
        ext = file.split('.')[-1]
    
    coord_class = import_module('.'+ext+'.'+ext.title())
    return coord_class(file, mode, buffer)
    