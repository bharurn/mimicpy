def readLine(line):
    vals = {}
    vals['record'] = line[:6].strip()
    
    if vals['record'] != 'HETATM' and vals['record'] != 'ATOM':
        vals['content'] = line[6:]
        return vals
    
    vals['serial'] = line[6:11].strip()
    vals['name'] = line[12:16].strip()
    vals['altLoc'] = line[16]
    vals['resName'] = line[17:20].strip()
    vals['chainID'] = line[21]
    vals['resSeq'] = line[22:26].strip()
    vals['iCode'] = line[26]
    vals['x'] = line[30:38].strip()
    vals['y'] = line[38:46].strip()
    vals['z'] = line[46:54].strip()
    vals['occupancy'] = line[54:60].strip()
    vals['tempFactor'] = line[60:66].strip()
    vals['element'] = line[76:78].strip()
    vals['charge'] = line[78:80].strip()
    
    return vals
