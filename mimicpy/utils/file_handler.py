from .errors import MiMiCPyError

def read(file, mode='r'):
    if mode not in ['r', 'rb']:
        raise MiMiCPyError('Mode {} is not a valid read mode. Only r and rb are allowed'.format(mode))
    with open(file, mode) as f:
        out = f.read()
    return out

def write(out, file, mode='w'):
    if mode not in ['w', 'wb']:
        raise MiMiCPyError('Mode {} is not a valid write mode. Only w and wb are allowed'.format(mode))
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
        self.close()

    def _next(self):
        if self.is_closed:
            raise StopIteration()
        out = self.file.read(self.buffer)
        if out == '':
            self.close()
            raise StopIteration()
        return out

    def close(self):
        self.is_closed = True
        if hasattr(self, 'file'):
            self.file.close()

    def readline(self):
        return self.file.readline()