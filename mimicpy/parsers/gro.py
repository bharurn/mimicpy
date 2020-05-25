from .._global import _Global as _global

def getBox(gro):
    # fast access of last line, requires UNIX shell
    tail = _global.host.run(f'tail -n 1 {gro}')
    return [float(v) for v in tail.split()]