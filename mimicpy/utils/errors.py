class MiMiCPyError(Exception):
   """Generic exception from MiMiCPy"""

class ExecutionError(MiMiCPyError):
    
    def __init__(self, cmd, msg):
        super().__init__(msg)
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
        super().__init__(msg)
        self.msg = msg
        
    def __str__(self):
        return f"{self.cmd} not set! Please set it in host.{self.msg}"

class ScriptError(MiMiCPyError):
     def __init__(self, cmd, msg):
        super().__init__(msg)
        self.msg = msg
        
     def __str__(self):
        return f"{self.msg} not a paramater of the script"

def defaultHook(cmd, out):
    if 'error' in out.lower(): raise ExecutionError(cmd.split(';')[-1], out)