import pickle
from .._global import _Global as gbl
from . import _mpt_writer

def write(inp, out):
    tail = gbl.host.run(f'tail -n 30 {inp}')
    mols = _mpt_writer.molecules(tail)
        
    file = gbl.host.open(inp, 'rb')
    atm_types_to_symb = _mpt_writer.atomtypes(file)
    ap = _mpt_writer.AtomsParser(file, mols, atm_types_to_symb)
    
    pkl = gbl.host.open(out, 'wb')
    pickle.dump(ap.mol_df, pkl)
    
class Reader:
    def __init__(self, file):
        pkl = gbl.host.open(file, 'rb') # open as bytes
        self.mpt = pickle.load(pkl) # unpickle
        pkl.close()
    
    def selectAtom(self, idx, mol=None):
        if mol:
            df = self.mpt[mol]
            return df[idx]
        else:
            for k,v in self.mpt.items():
                if idx > v[2]: continue
                else:
                    return df[idx]