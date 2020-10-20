from .errors import MiMiCPyError

def read(file, mode):
    if mode not in ['r', 'rb']:
        raise MiMiCPyError(f'Mode {mode} is not a valid read mode. Only r and rb are allowed.')
    with open(file, mode) as f:
        out = f.read()
    return out

def write(out, file, mode):
    if mode not in ['w', 'wb']:
        raise MiMiCPyError(f'Mode {mode} is not a valid write mode. Only w and wb are allowed.')
    with open(file, mode) as f:
        f.write(out)
        
class Parser:
    """implements methods for iterable objects and wraps around readline"""

    def __init__(self, file, buffer=1000):
        self.file = open(file, 'r')
        self.buffer = buffer
        self.is_closed = False

    def __iter__(self):
        return self

    def __next__(self):
        return self._next()

    def __del__(self):
        self._del()

    def _next(self):
        if self.is_closed:
            raise StopIteration()
        out = self.file.read(self.buffer)
        if out == '':
            self._del()
            raise StopIteration()
        return out

    def _del(self):
        self.is_closed = True
        if hasattr(self, 'file'):
            self.file.close()

    def readline(self):
        return self.file.readline()