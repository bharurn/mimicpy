import pytest
import mimicpy
import pickle

def test_dppc():
    mimicpy.setHost('dppc')
    mol_list, topol_dict = mimicpy.parsers.top.read('topol.top', guess_elems=False)
    assert mol_list==[('DPPC', 2), ('SOL', 100)], 'mol_list no parsed correctly!'
    with open('result.pkl', 'rb') as f:
        known_topol_dict = pickle.load(f)
    assert topol_dict.repeating == known_topol_dict.repeating
    assert topol_dict.dict_df.keys() == known_topol_dict.dict_df.keys()
    for k,v in topol_dict.dict_df.items():
        assert v.equals(v)