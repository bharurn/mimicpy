from . import top_reader
import re
import os

def read(topol_file, nonstd_atm_types={}, buff=1000, guess_elems=True):
    file = top_reader.Parser(topol_file, buff)
    topol_txt = "".join(file)
    include_file_list = re.compile(r"#include\s+[\"\'](.+)\s*[\"\']", re.MULTILINE).findall(topol_txt)
                                   
    ffnonbonded = os.path.join(os.path.dirname(include_file_list[0]), 'ffnonbonded.itp')
    atm_types = top_reader.atomtypes(ffnonbonded, buff)
    atm_types.update(nonstd_atm_types)
    
    mols_data = top_reader.molecules(topol_txt)
    mols = [m[0] for m in mols_data]
    
    top_reader.ITPParser.mol_df = {}
    
    itp_parser = top_reader.ITPParser(mols, atm_types, buff, guess_elems)
    
    itp_parser.parse('topol.top', topol_txt)
    
    for file in include_file_list[1:]:
        itp_parser.parse(file)
    
    mol_df = dict(zip(itp_parser.mols, itp_parser.dfs))
    
    return mols_data, mol_df