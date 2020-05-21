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
        writer = getattr(self, option)
        if writer != None:
            writer.write(value+'\n')
            writer.flush()
    
    def __del__(self):
        for k in self.kwargs:
            v = getattr(self, k)
            v.close()
                
    def close(self):
        self.__del__()
        
class LogString:
    def __init__(self):
        self.val = ''
    
    def write(self, s):
        self.val += s
    
    def flush(self):
        pass
    
    def close():
        pass
    
    def __str__(self):
        return self.val
    
class StdOut:
    def write(self, s):
        sys.stdout.write(s)
    
    def flush(self):
        pass
    
    def close(self):
        pass