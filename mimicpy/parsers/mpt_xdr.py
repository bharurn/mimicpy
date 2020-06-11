from .top import TopolDict
from .top_reader import ITPParser
import pandas as pd

encoder = 'utf-8'

def pack_strlist(packer, lst):
    """Pack list of strings as comma seperated string"""
    s = ",".join(lst).encode(encoder)
    packer.pack_string(s)

def pack_df(packer, df):
    # pack index first
    packer.pack_list(df.index.to_list(), packer.pack_int)
    for i, col in enumerate(df.columns): #iterate thru columns
        lst = df[col].to_list() #convert column to list
        if i in [0,2,3,5]: #pack string columns, .i.e, atom name, type, resname, etc.
            pack_strlist(packer, lst)
        elif i == 1: #pack resid column as list of ints
            packer.pack_list(lst, packer.pack_int)
        elif i in [4,6]: #pack charge and mass as list of floats
            packer.pack_list(lst, packer.pack_float)

def pack_topol_dict(packer, topol_dict):
    # pack repearting dict
    pack_strlist(packer, topol_dict.repeating.keys()) #keys
    pack_strlist(packer, topol_dict.repeating.values()) #values
    
    # pack dict of dataframes
    for k,v in topol_dict.dict_df.items():
        packer.pack_string(k.encode('utf-8')) #dict key/name of mol
        pack_df(packer, v) #dataframe/topology


def unpack_strlist(unpacker):
    return unpacker.unpack_string().decode(encoder).split(',')

def unpack_df(unpacker):
    lst1 = unpacker.unpack_list(unpacker.unpack_int) # unpack atom number list
    lst2 = unpack_strlist(unpacker) # unpack atom types list
    lst3 = unpacker.unpack_list(unpacker.unpack_int) # resid list
    lst4 = unpack_strlist(unpacker) # unpack resnames
    lst5 = unpack_strlist(unpacker) # unpack atom names
    lst6 = unpacker.unpack_list(unpacker.unpack_float) # unpack charges
    lst7 = unpack_strlist(unpacker) # unpack elemtn
    lst8 = unpacker.unpack_list(unpacker.unpack_float) # unpack masses
    df = pd.DataFrame([lst1, lst2, lst3, lst4, lst5, lst6, lst7, lst8]).T
    df.columns = ITPParser.columns
    return df.set_index(df.columns[0])

def unpack_topol_dict(unpacker):
    # unpack repeating dict
    repeating_keys = unpack_strlist(unpacker)
    repeating_vals = unpack_strlist(unpacker)
    repeating = dict(zip(repeating_keys, repeating_vals))
    
    dict_df = {}
    while True: # unpack dataframe dict until EOF
        try:
            mol = unpacker.unpack_string().decode(encoder)
        except EOFError:
            break
        df = unpack_df(unpacker)
        dict_df[mol] = df
    
    return TopolDict(dict_df, repeating)