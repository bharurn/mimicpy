from .._global import _Global as gbl
import numpy as np
import pandas as pd

def getBox(gro):
    # fast access of last line, requires UNIX shell
    tail = gbl.host.run(f'tail -n 1 {gro}')
    return [float(v) for v in tail.split()]

def read(file, lines=500):
    d = gbl.host.open(file, 'rb')
    
    d.readline()
    no = int(d.readline().decode())
    
    s = d.readline().decode()
    buff = len(s)*lines
    i = 0
    
    def mapper(a):
        if a.isnumeric(): return np.nan
        try:
            return float(a)
        except:
            return np.nan
    
    vals = np.array(list(map(mapper, s.split())))
    
    rows = len(vals[~np.isnan(vals)]) # 6--> force, 3--> only coords
    
    while i < no:
        s = d.read(buff).decode()
        i += lines
        
        vals = np.append(vals, np.array(list(map(mapper, s.split()))) )
    
    vals = vals[~np.isnan(vals)]
    
    expected_len = no*rows
    
    # checking if box size is read already
    if len(vals) == expected_len:
        # box vector not been read yet
        coords = vals
        box = d.readline().decode()
    elif len(vals) == expected_len+3:
        # box vector already read
        coords = vals[:-3]
        box = vals[-3:]
    else:
        #raise ParserError
        pass
    
    # check whether forces present or no before reshaping
    if rows == 6: cols = ['x', 'y', 'z', 'force-x', 'force-y', 'force-z']
    elif rows == 3: cols = ['x', 'y', 'z']
    
    coords = pd.DataFrame(coords.reshape(rows, no).T, columns=cols)
    coords['id'] = coords.index.to_numpy()+1
    
    return coords, box.tolist()