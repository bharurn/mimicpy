from ..__init__ import host
import re
import pandas as pd
from .plotbox import plotBoxGen

def xvg(cls, xvg):
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
    
    df = pd.DataFrame([x,y], columns=(xlabel, ylabel))
    
    return df.set_index([xlabel])
    
def gmxErrors(file_eval=None):
    
    def _parse(out):
        output = host.cmd.read(out)   
        return host._notes(output), host._errorHandle(output, dont_raise=True)
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

@plotBoxGen('Logs', 'Step', 'Time')
def log(file):
    f = host.cmd.read(file)
    x = re.compile(r"^\s*(Step\s*Time)\n(.+)\n(?:\n.+)+\s*Energies\s*\(kJ/mol\)((?:\n.+)+)\s+\n",\
               re.MULTILINE)
    res = x.findall(f)
    
    def c(log_txt):
        colms = []
        vals = []
        n = 15
        colms += log_txt[0].split()
        vals += [float(i) for i in log_txt[1].split()]
        for j,l in enumerate(log_txt[2].splitlines()):
            if j%2 != 0:
                cols = [l[i:i+n].strip() for i in range(0, len(l), n)]
                colms.extend(cols)
            else:
                cols = [float(l[i:i+n].strip()) for i in range(0, len(l), n)]
                vals.extend(cols)
        return (colms, vals)
   
    cols, v1 = c(res[0])
    
    vals = []
    for i in res:
        vals.append(c(i)[1])
        
    df = pd.DataFrame(vals, columns=cols)
    df = df.set_index(['Step', 'Time'])
    
    return df