"""Module for custom exceptions"""


class MiMiCPyError(Exception):
    """Generic exception from MiMiCPy"""

class SelectionError(MiMiCPyError):
    """Error in selecting atoms in Mpt"""

class ParserError(MiMiCPyError):
    """Error in parsing a file"""
    def __init__(self, file='', file_type='', line_number='', details=''):
        self.file = file
        self.file_type = file_type
        self.line_number = line_number
        self.details = details
        if self.file_type:
            self.file_type = ' as ' + file_type
        if self.line_number:
            self.line_number = f' at line number {line_number}'
        if self.details:
            self.details = ': ' + details

    def __str__(self):
        return f'Error parsing {self.file}{self.file_type}{self.line_number}{self.details}.'

class ScriptError(MiMiCPyError):
    """Requested parameter does not exist in script object"""
    def __init__(self, parameter):
        self.parameter = parameter

    def __str__(self):
        return f'The {self.parameter} parameter has not been set or has been set incorrectly.'

def asserter(boolean, error, *args, **kwargs):
    """Custom assert to raise custom errors"""
    if not boolean:
        raise error(*args, **kwargs)
