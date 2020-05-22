from .._global import _Global as _global

def getBox(gro):
    last = _global.host.run(f'tail -n 1 {gro}') # fast access of last line, requires UNIX shell
    
    return [float(v) for v in last.split()]