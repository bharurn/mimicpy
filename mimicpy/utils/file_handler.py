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