from mimicpy.parsers.mpt import MPT

def test_translate():
    mpt = MPT([], {}, 'w') # empty MPT object
    ev, keys = mpt._MPT__translate('resname is SER or resid is 132')
    assert ev == "(np_vals['resname']=='SER') | (np_vals['resid']==132)"
    assert keys == ['resname', 'resid']
    
    ev, keys = mpt._MPT__translate('resname is HIS or resid is 100 and name is CA')
    assert ev == "(np_vals['resname']=='HIS') | (np_vals['resid']==100) & (np_vals['name']=='CA')"
    assert keys == ['resname', 'resid', 'name']

def getFakeMPT():
    mpt = MPT([], {}, 'w') # empty MPT object
    mpt._MPT__mot_list = [('MOL1', 1)]
    mpt._MPT__all_data = [['CA', 'CA', 'CT', 'H', 'N*'],\
                          [1, 1, 1, 2, 2], \
                          ['UNK1', 'UNK1', 'UNK1', 'UNK2', 'UNK2'], \
                          ['C1', 'C2', 'C3', 'H1', 'N1'], \
                          [0.5, -0.5, 0.2, 0, 2], \
                          ['C', 'C', 'C', 'H', 'N'], \
                          [12, 12, 12, 1, 14]]
    return mpt

def test_prop():
    mpt = getFakeMPT()
    
    assert mpt['type'] == ['CA', 'CA', 'CT', 'H', 'N*']