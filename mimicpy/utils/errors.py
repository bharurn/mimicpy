# This is part of MiMiCPy

"""

This module contains classes for all custom exceptions
used in the package

"""

class MiMiCPyError(Exception):
   """Generic exception from MiMiCPy"""

class ParserError(MiMiCPyError):
    """
    Raised when a file could not be parsed
    using the format requsted
    """
    def __init__(self, ftype, file=None, no=None):
        if file is None: file = 'data block'
        else: file = 'file ' + file
        self.file = file
        self.ftype = ftype
        self.no = no
        
    def __str__(self):
        s = f"{self.file} cannot be parsed as {self.ftype}"
        if self.no:
            s = f"Line number {self.no} in "+s
        else:
            s = s[0].upper() + s[1:] # make data/file --> Data/File if it is first letter
            
        return s

class ExecutionError(MiMiCPyError):
    """
    Raised when a command executed by host exits
    with and error message
    """
    def __init__(self, cmd, msg):
        self.cmd = cmd
        self.msg = msg
        
    def __str__(self):
        return f"Command attempted {self.cmd}\n{self.msg}"

class GromacsError(ExecutionError):
    """
    Raised when a Gromacs command executed by host
    exits with and error message
    """

class SlurmBatchError(ExecutionError):
    """
    Raised when a Slurm sbatch command executed by host
    exits with and error message
    """
    def __str__(self):
        return f"Attempted to run sbatch {self.cmd}\nself.msg"

class EnvNotSetError(MiMiCPyError):
    """
    Raised when an enviornment path requested by
    MiMiCPy is not set
    """
    def __init__(self, cmd, msg):
        self.msg = msg
        
    def __str__(self):
        return f"{self.cmd} not set! Please set it in host.{self.msg}"

class ScriptError(MiMiCPyError):
    """
    Raised when a requested parameter of a script
    object does not exist
    """
    def __init__(self, cmd, msg):
        self.msg = msg
        
    def __str__(self):
        return f"{self.msg} not a paramater of the script"

def defaultHook(cmd, out):
    """Default error checking hook called by host"""
    if 'error' in out.lower() or 'not recognisable' in out.lower(): raise ExecutionError(cmd.split(';')[-1], out)
    
def asserter(boolean, error, *args, **kwargs):
    """Custom assert to raise custom errors"""
    if not boolean: raise error(*args, **kwargs)