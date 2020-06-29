from . import top_reader
from .._global import _Global as gbl
from . import pdb as parse_pdb

class TopolDict:
    def __init__(self, dict_df, repeating):
        self.dict_df = dict_df
        self.repeating = repeating
    
    @classmethod
    def fromDict(cls, df):
        keys = list(df.keys())
        df2 = df.copy()
        repeating = {}
        for i in range(len(keys)):
            key_i = keys[i]
            for j in range(i+1, len(keys)):
                key_j = keys[j]
                if df[key_i].equals(df[key_j]):
                    repeating[key_j] = key_i
                    del df2[key_j]
        return cls(df2, repeating)
    
    def __getitem__(self, key):
        if key in self.dict_df:
            return self.dict_df[key]
        elif key in self.repeating:
            return self.dict_df[self.repeating[key]]
        else:
            raise KeyError(f"Molecule {key} not in topology")
    
    def __getAll(self):
        extras = self.dict_df.copy()
        for i in self.repeating:
            extras[i] = self.__getitem__(i)
        return extras
    
    def __str__(self): return str(self.__getAll())

    def __repr__(self): return repr(self.__getAll())

def read(topol_file, nonstd_atm_types={}, buff=1000, guess_elems=True):
    file = top_reader.Parser(topol_file, buff)
    # we assume that .top is not too big and can be fully loaded to memory
    # gromacs usually only writes at most one molecule to .top file
    topol_txt = "".join(file)
    
    # get all itp files in .top
    include_file_list = top_reader.include_file_regex.findall(topol_txt)
    
    dirname = gbl.host.dirname(topol_file) # get dir of topol file
    include_file_list = [gbl.host.join(dirname, i) for i in include_file_list]
    
    # look for atomtypes in first itp
    atm_types = top_reader.atomtypes(include_file_list[0], buff)
    atm_types.update(nonstd_atm_types)
    
    mols_data = top_reader.molecules(topol_txt) # mol, no. list
    mols = [m[0] for m in mols_data]
    
    # clear the mol and dfs of ITPParser
    # these are static vars of class, so is shared b/w objects
    top_reader.ITPParser.clear()
    
    itp_parser = top_reader.ITPParser(mols, atm_types, buff, guess_elems) # init itp parser
    
    # parse .top file, in case some atoms defined there
    itp_parser.parse('topol.top', topol_txt)
    
    # parse itp files
    for file in include_file_list[1:]:
        if gbl.host.fileExists(file):
            itp_parser.parse(file)
        else:
            gbl.logger.write('warning', f"Cannot find {file}, skipping..")
    
    mol_df = dict(zip(itp_parser.mols, itp_parser.dfs))
    return mols_data, TopolDict.fromDict(mol_df)

def nonStdTypes(pdb, itp, *resnames, buff=1000):
    pdb_df = parse_pdb.parseFile(pdb, buff)
    
    # get chain info
    chains = pdb_df['chainID'].unique().tolist()
    chains.remove(' ') # remove elements wtih no chain info from lisy
    
    # either select only first chain of residue or all residues with no chain info
    pdb_df = pdb_df[(pdb_df['chainID'] == chains[0]) | (pdb_df['chainID'] == ' ')]
    
    # get elems
    elems = [a for name in resnames for a in pdb_df['element'][pdb_df['resname']==name].to_list()]
    
    # get atom types
    top_reader.ITPParser.clear()
    
    # fake class to satisfy atomtypes dict
    class defdict:
        def __contains__(self, item): return True
        def __getitem__(self,key): return ' '
    
    itp_parser = top_reader.ITPParser(resnames, defdict(), buff, False)
    itp_parser.parse(itp)
    atm_types = [a for df in itp_parser.dfs for a in df['type'].to_list()]
    
    # assert len(atm_types) == len(elems) before zipping
    return dict(zip(atm_types, elems))