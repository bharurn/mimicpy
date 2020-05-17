from collections import OrderedDict 
from ..utils.errors import ScriptError

class Script(object):
    def __init__(self, **kwargs):
        super().__setattr__('__orddict__', OrderedDict())
        
        for key, value in kwargs.items():
            setattr(self, key, value)
            
    def __setattr__(self, key, value):
        if key.startswith('_') or key == 'hasparams' or key == 'params':
            self.__dict__[key] = value
        else:
            self.__orddict__[key] = value
        
    def __getattr__(self, key):
        if key.startswith('_') or key == 'hasparams' or key == 'params':
            return self.__getattribute__('__dict__')[key]
        else:
            if key not in self.__getattribute__('__orddict__'):
                raise ScriptError(key)
                
            return self.__getattribute__('__orddict__')[key]
    
    def hasparam(self, key):
        if key in self.__orddict__:
            return True
        else:
            return False
    
    def params(self): return self.__orddict__