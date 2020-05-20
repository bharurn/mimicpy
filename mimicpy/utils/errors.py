class MiMiCPyError(Exception):
   """Generic exception from MiMiCPy"""

class ParserError(MiMiCPyError):
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
    
    def __init__(self, cmd, msg):
        self.cmd = cmd
        self.msg = msg
        
    def __str__(self):
        return f"Command attempted {self.cmd}\n{self.msg}"

class GromacsError(ExecutionError):
    pass

class SlurmBatchError(ExecutionError):
    def __str__(self):
        return f"Attempted to run sbatch {self.cmd}\nself.msg"

class EnvNotSetError(MiMiCPyError):
    def __init__(self, cmd, msg):
        self.msg = msg
        
    def __str__(self):
        return f"{self.cmd} not set! Please set it in host.{self.msg}"

class ScriptError(MiMiCPyError):
     def __init__(self, cmd, msg):
        self.msg = msg
        
     def __str__(self):
        return f"{self.msg} not a paramater of the script"

def defaultHook(cmd, out):
    if 'error' in out.lower(): raise ExecutionError(cmd.split(';')[-1], out)
    
def asserter(boolean, error, *args, **kwargs):
    if not boolean: raise error(*args, **kwargs)