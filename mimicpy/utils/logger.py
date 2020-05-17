import sys
from .errors import MiMiCPyError

class Logger:
    def __init__(self, **kwargs):
       
        self.kwargs = kwargs
        
        for k,v in kwargs.items():
            setattr(self, k, v)
        
    def write(self, option, value):
        if option not in self.kwargs:
            raise MiMiCPyError(f'{option} is not a logger stream')
        writer = self._get(option)
        if writer != None:
            writer.write(value+'\n')
    
    def __del__(self):
        for k in self.kwargs:
            v = self._get(k)
            if v != None and v != sys.stdout:
                v.close()
                
    def close(self):
        self.__del__()