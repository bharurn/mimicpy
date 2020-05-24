# This is part of MiMiCPy

"""

This module contains the Logger, LogString and StdOut
classes to allow for logging multiple streams and redirection

"""

import sys
from .errors import MiMiCPyError

class Logger:
    def __init__(self, **kwargs):
       
        self.kwargs = list(kwargs.keys())
        
        for k,v in kwargs.items():
            setattr(self, k, v)
    
    def add(self, **kwargs):
       
        self.kwargs.extend(list(kwargs.keys()))
        
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
    
    def close(self):
        pass
    
    def __str__(self):
        return self.val
    
    def __repr__(self):
        return self.val
    
class StdOut:
    def write(self, s):
        print(s, end='')
    
    def flush(self):
        pass
    
    def close(self):
        pass