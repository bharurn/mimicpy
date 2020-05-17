class MiMiCPyError(Exception):
   """Base class for other exceptions"""
   pass

class RemoteShellError(MiMiCPyError):
    pass

class OpenBabelError(MiMiCPyError):
    pass

class LEaPError(MiMiCPyError):
    pass

class AcpypeError(MiMiCPyError):
    pass

class GromacsError(MiMiCPyError):
    pass

class SlurmError(MiMiCPyError):
    pass

class EnvNotSetError(MiMiCPyError):
    pass

class ScriptError(MiMiCPyError):
    pass