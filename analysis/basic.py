from pygmx import host
from pygmx.run.gmxrun import GMX
import re

class Value:
    def __init__(self, xlabel, ylabel, x, y): self.xlabel=xlabel; self.ylabel=ylabel; self.x=x; self.y=y
    
    def __add__(self, a):
        xlabel = a.xlabel
        ylabel = a.ylabel
        x = self.x  + a.x
        y = self.y + a.y
        
        return Value(xlabel, ylabel, x, y)
    
    @classmethod
    def empty(cls):
        return cls('', '', [], [])
    
    @classmethod
    def fromxvg(cls, xvg):
        xlabel = 'X Axis'
        ylabel = 'Y Axis'
        x = []
        y = []
    
        reg = re.compile(r"@\s*xaxis\s*label\s*\"(.*?)\"", re.MULTILINE).findall(xvg)
        if reg == []: xlabel = 'X Axis'
        else: xlabel = reg[0]
    
        reg = re.compile(r"@\s*yaxis\s*label\s*\"(.*?)\"", re.MULTILINE).findall(xvg)
        if reg == []: ylabel = 'Y Axis'
        else: ylabel = reg[0]
    
        reg = re.compile(r"@\s*s0\s*legend\s*\"(.*?)\"", re.MULTILINE).findall(xvg)
        if reg != []: ylabel = reg[0] + ' ' + ylabel
    
        for line in xvg.splitlines():
            if line.startswith('#') or line.startswith('@'):
                pass
            else:
                splt = line.split()
                x.append(float(splt[0]))
                y.append(float(splt[1]))
                
        return cls(xlabel, ylabel, x, y)

def _parse(out):
        
    output = host.cmd.read(out)   
    return host._notes(output), host._errorHandle(output, dont_raise=True)
    
def parseOutput(file_eval=None):
    
    # include ls() in local
    if file_eval is None:
        files = host.cmd.ls(file_eval=lambda a: True if a.endswith('.log') or a.endswith('.out') else False)
    else:
        files = host.cmd.ls(file_eval=file_eval)
        
    for file in files:
        print(f"\n<========{file}========>")
            
        notes, errors = _parse(file)
            
        if notes: print(f"\nThere were note(s)/warning(s):\n\n{notes}")
        if errors: print(f"\nThere were error(s):\n\n{errors}")

def dump(file):
    ext = file.split('.')[-1]
    if ext == 's': cmd = 's'
    elif ext == 'trr' or ext == 'xtc' or ext == 'gro' or ext == 'pdb': cmd = 'f'
    elif ext == 'edr': cmd = 'e'
    elif ext == 'cpt': cmd = 'cp'
    elif ext == 'top': cmd = 'p'
    elif ext == 'mtx': cmd = 'mtx'
    
    out = host.cmd.run(f'{host.gmx} dump -{cmd} {file}')
    return out

def energy(file, component, **kwargs):
    
    txt = GMX.gmx('energy', f = file, **kwargs, noverbose=True, stdin=component)
    
    try:
        out = host.cmd.read('energy.xvg')
    except:
        raise Exception(f'Error {txt}')
        
    host.cmd.rm('energy.xvg')
    
    return Value.fromxvg(out)