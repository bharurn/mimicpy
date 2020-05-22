from collections import defaultdict, OrderedDict
import pandas as pd
from ..utils.constants import element_names
 
def isSection(section, txt):
    if section == '*': section = ''
    if '[' in txt and ']' in txt and section in txt:
        return True
    else: return False
    
def molecules(tail):    
    _mols = defaultdict(int)
        
    for line in tail.splitlines()[::-1]:
        if isSection('molecules', line):
            break
        elif line.strip() == '' or line.startswith(';'):
            continue
        mol, no = line.split()
        _mols[mol] += int(no)
            
    return OrderedDict(list(_mols.items())[::-1])
    
def atomtypes(f, buff):
    atomtypes = ''
    while True:
        s = f.read(buff).decode()
        atomtypes += s
        if isSection('bondtypes', s): break
        
    atm_types_to_symb = {} # init atom types to symbol
    start = False
    
    for line in atomtypes.splitlines():
        if isSection('atomtypes', line) and atm_types_to_symb == {}: # read only first atomtypes
            start = True
        elif start:
            if isSection('*', line):
                start = False
            elif not line.startswith(';') and line.strip() != '': # ignore comments and null lines
                e = line.split()[0] # first val is atom type
                _n = line.split()[1]
                
                if not _n.isnumeric(): continue # check just in case
                else: n = int(_n)-1 # second is element no.
                
                if n == -1: continue # dummy masses, skip for now
                atm_types_to_symb[e]  = element_names[n] # fill up at
    
    return atm_types_to_symb


class AtomsParser:

    def parseAtoms(self, buff):
        start = False
        atoms = ''
        end = ['bonds']
        while True:
            s = self.f.read(buff).decode()
            if isSection('moleculetype', s):
                start = True
            if start:    
                atoms += s
                if any([isSection(hdr, s) for hdr in end]) or s.count('moleculetype') > 1:
                    start = False
                    break
                if s == '':
                    break
        
        df_ = {'number':[], 'type':[], 'resNo':[], 'resName':[], 'name':[], 'charge':[], 'element':[], 'mass':[]}
        start = False
        mol = ''
        end.append('moleculetype')
        atom_splt = atoms.splitlines()
        for i, line in enumerate(atom_splt):
            if any([isSection(hdr, line) for hdr in end]):
                start = False
                if mol in self.mol_df.keys():
                    self.mol_df[mol][2] = self.natms
                    self.mol_df[mol][3] = pd.DataFrame(df_).set_index(['number'])
                    self.mol_df[mol][1] = len(self.mol_df[mol][3])
                    self.natms += self.mol_df[mol][0]*self.mol_df[mol][1]
                    df_ = {k:[] for k in df_.keys()}
                if isSection('moleculetype', line):
                    for l in atom_splt[i+1:]:
                        if not l.startswith(';'):
                            mol = l.split()[0]
                            print(mol)
                            break
            if isSection('atoms', line) and mol in self.mol_df.keys(): # read all atom sections
                start = True
            elif start and not line.startswith(';') and line.strip() != '': # ignore comments and null lines
                splt = line.split()
                if len(splt) >= 8:
                    nr, _type, resnr, res, name, cgnr, q, mass = splt[:8]
                else:
                    continue
                df_['number'].append(int(nr))
                df_['type'].append(_type)
                df_['resNo'].append(int(resnr))
                df_['resName'].append(res)
                df_['name'].append(name)
                df_['charge'].append(float(q))
                
                if _type in self.atm_types_to_symb:
                    elem = self.atm_types_to_symb[_type]
                else:
                    elem = ''
                
                df_['element'].append(elem)
                df_['mass'].append(mass)
    
    def __init__(self, f, mols, atm_types_to_symb, buff):
        self.f = f
        self.atm_types_to_symb = atm_types_to_symb
        # mol_name: (no. of mols, no. of atoms in one mol, no. of atoms before mol, df)
        self.mol_df = OrderedDict( (k,[v,0,0,pd.DataFrame()]) for k,v in mols.items())
        self.natms = 0
        
        while any(m[3].empty for k,m in self.mol_df.items()): 
            self.parseAtoms(buff)
