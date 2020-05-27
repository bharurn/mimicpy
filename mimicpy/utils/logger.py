# This is part of MiMiCPy

"""

This module contains the Logger, LogString and StdOut
classes to allow for logging multiple streams and redirection

"""
from .._global import _Global as gbl
from .errors import MiMiCPyError

class Logger:
    def __init__(self, **kwargs):
       
        self.kwargs = list(kwargs.keys())
        
    def add(self, **kwargs):
       
        self.kwargs.extend(list(kwargs.keys()))
        
    def write(self, option, value):
        if option not in self.kwargs:
            raise MiMiCPyError(f'{option} is not a logger stream')
        writer = getattr(self, option)
        if writer != None:
            writer.write(value+'\n')
        
class LogString:
    def __init__(self):
        self.val = ''
    
    def write(self, s):
        self.val += s
    
    def __str__(self):
        return self.val
    
    def __repr__(self):
        return self.val
    
class LogFile:
    def __init__(self, name, forceLocal=False):
        self.fname = name
        self.forceLocal = forceLocal
        
    def write(self, s):
        if not self.forceLocal:
            gbl.host.write(s, self.fname)
        else:
            with open(self.fname, 'r') as f:
                f.write(s)