import pickle
from .._global import _Global as gbl
from . import _mpt_writer
import pandas as pd

def write(inp, out, nonstd_atm_types={}, buff=1000):
    tail = gbl.host.run(f'tail -n 30 {inp}')
    mols = _mpt_writer.molecules(tail)
        
    file = gbl.host.open(inp, 'rb')
    atm_types_to_symb = _mpt_writer.atomtypes(file, buff)
    # extend atm_types_to_symb with nonstdligands
    # should be dict of atom type --> symb
    atm_types_to_symb.update(nonstd_atm_types)
    ap = _mpt_writer.AtomsParser(file, mols, atm_types_to_symb, buff)
    
    df = ap.mol_df
    
    # replace repeating dataframes with the string name of prev mol
    keys = list(df.keys())
    vals = []
    for i in range(len(keys)):
        key_i = keys[i]
        for j in range(i+1, len(keys)):
            key_j = keys[j]
        
            try:
                if all(df[key_i][3] == df[key_j][3]): 
                    vals.append((key_i, key_j))
            except ValueError:
                continue
    
    for k,v in vals: df[v][3] = k
    
    pkl = gbl.host.open(out, 'wb')
    pickle.dump(df, pkl)
    
class Reader:
    
    def __init__(self, file):
        pkl = gbl.host.open(file, 'rb') # open as bytes
        self.mpt = pickle.load(pkl) # unpickle
        pkl.close()
    
    def _get_df(self, mol):
        third_val = self.mpt[mol][3]
        
        if isinstance(third_val, str):
            return self.mpt[third_val][3]
        else:
            return third_val
    
    def selectAtom(self, idx, mol=None, relative=False):
        orig_id = idx
        if not mol:
            mol = ''
            for k,v in self.mpt.items():
                if idx <= v[2]: break
                else: mol = k 
        
        df = self._get_df(mol)
        if not relative:
            atms_before = self.mpt[mol][2]
            idx = idx-atms_before
        
        natms = self.mpt[mol][1]
        
        if natms > 1:
            idx -= 1
            col = idx % natms
            idx = col+1
            
        srs = df.loc[idx]
        
        return srs.append(pd.Series({'mol':mol, 'id':orig_id}))
    
    def selectAtoms(self, ids):
        s = [self.selectAtom(i) for i in ids]
        return pd.concat(s, axis=1).T
    
    # TO DO: make func to merge mpt with gro file
    
    # TO DO: need to rewrite this for current mpt format
    def parse_selec(selection, df):
        """Translate selection language string into pandas dataframe selection"""
    
        # if selection is a lambda func, then just call and return it
        # this is provided for debuggin puposes
        LAMBDA = lambda:0
        if isinstance(selection, type(LAMBDA)):
            return df[selection(df)]
            
        # below code translates selection langauge to pandas boolean
        # selection eg., resName is SER and number < 25 and chainID not B
        # will be translated to df['resName'] == 'SER' and df.index == 25 and df['chainID'] != 'B'
    
        ev = '' # converted string
        i = 0 # counter to keep track of word position
        for s in selection.split():
            if i == 0: # if starting of set
                ev += f"(df['{s}']"
                # if and/or encountered, reset i to -1 (will become 0 due to i+= 1 at end)
                # so we can start parsing again
            elif s == 'or':
                ev += f' | '
                i = -1
            elif s == 'and':
                ev += f' & '
                i = -1
            elif s == 'is':
                ev += '=='
            elif s == 'not':
                ev += '!='
            elif s == '>' or s == '>=' or s == '<' or s == '<=':
                ev += s
            else: # parse everything else, meant for the third word
                if s.isnumeric():
                    ev += f"{s})"
                else:
                    ev += f"'{s}')"
                
            i += 1

        ev = f"df.loc[{ev}]" # eg., df.loc[ df['resName'] == 'SER' and df.index == 25 and df['chainID'] != 'B' ]
        ev = ev.replace("df['number']","df.index") # replace df['number'] to df.index as number is the index of df
        gbl.logger.write('debug2', f'Selection command translated to: {ev}')
        
        return eval(ev) # evaluate string and return the dataframe