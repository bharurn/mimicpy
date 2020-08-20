""" Module for custom exceptions """


class MiMiCPyError(Exception):
    """ Generic exception from MiMiCPy """


class SelectionError(MiMiCPyError):
    """ Error in selecting atoms in Mpt """


class ParserError(MiMiCPyError):
    """ Error in paring a file """
    def __init__(self, file='', file_type='', line_number='', details=''):
        self.file = file
        self.file_type = file_type
        self.line_number = line_number
        self.details = details

        if self.file_type:
            self.file_type = ' as ' + file_type
        if self.line_number:
            self.line_number = f" at line number {line_number}"
        if self.details:
            self.details = ': ' + details

    def __str__(self):
        error_message = f"Error parsing {self.file}{self.file_type}{self.line_number}{self.details}"
        return error_message

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

class EnvNotSetError(MiMiCPyError):
    """
    Raised when an enviornment path requested by
    MiMiCPy is not set
    """
    def __init__(self, env, cmd):
        self.env = env
        self.cmd = cmd

    def __str__(self):
        return f"{self.env} not set! Please set it with the keyword {self.cmd}."

class ScriptError(MiMiCPyError):
    """
    Raised when a requested parameter of a script
    object does not exist
    """
    def __init__(self, msg):
        self.msg = msg

    def __str__(self):
        return f"{self.msg} not a paramater of the script"

def defaultHook(cmd, out):
    """Default error checking hook called by host"""
    if 'error' in out.lower() or 'not recognisable' in out.lower() or 'not found' in out.lower():
        raise ExecutionError(cmd.split(';')[-1], out)

def asserter(boolean, error, *args, **kwargs):
    """Custom assert to raise custom errors"""
    if not boolean:
        raise error(*args, **kwargs)
