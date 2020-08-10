from mimicpy.io.mpt import Mpt
import pytest

def getMockTopol():
    import pandas as pd
    
    df1 = pd.DataFrame({'type': ['CA', 'CA', 'CT', 'H', 'N*'],\
                       'resid': [1, 1, 1, 2, 2],\
                       'resname': ['UNK1', 'UNK1', 'UNK1', 'UNK2', 'UNK2'],\
                       'name': ['C1', 'C2', 'C3', 'H1', 'N1'],\
                       'charge': [0.5, -0.5, 0.2, 0, 2],\
                       'element': ['C', 'C', 'C', 'H', 'N'],\
                       'mass': [12, 12, 12, 1, 14]})
    df1['number'] = df1.index+1
    
    df2 = pd.DataFrame({'type': ['NA'],\
                       'resid': [1],\
                       'resname': ['NA+'],\
                       'name': ['NA'],\
                       'charge': [1],\
                       'element': ['Na'],\
                       'mass': [23]})
    
    df2['number'] = df2.index+1
    
    return df1, df2
   
def test_topol_dict():
    from mimicpy.io.topol_dict import TopolDict
    df1, df2 = getMockTopol()
    topol_dict = TopolDict.from_dict({'MOL1':df1, 'MOL2':df1, 'NA1': df2, 'NA2':df2})
    
    assert topol_dict.repeating == {'MOL2': 'MOL1', 'NA2': 'NA1'}
    
    topol_dict = TopolDict.from_dict({'MOL1':df1, 'MOL2':df1, 'NA1': df2, 'NA2':df2, 'MOL3':df1, 'NA3': df2, 'NA4': df2})
    
    assert topol_dict.repeating == {'MOL2': 'MOL1', 'MOL3': 'MOL1', 'NA2': 'NA1', 'NA3': 'NA1', 'NA4': 'NA1'}

def test_prop():
    df1, df2 = getMockTopol()
    
    from mimicpy.io.topol_dict import TopolDict
    topol_dict = TopolDict.from_dict({'MOL1':df1, 'NA1':df2, 'NA2':df2, 'MOL2':df1})
    
    mpt = Mpt([('MOL1', 2), ('NA1', 1), ('MOL2', 1), ('NA2', 3), ('NA1', 1)], topol_dict, 'r')
    
    assert len(mpt['type']) == 20
    assert mpt['type'] == ['CA', 'CA', 'CT', 'H', 'N*']*2 + ['NA'] + ['CA', 'CA', 'CT', 'H', 'N*'] + ['NA']*4
    assert mpt['resname'] == (['UNK1']*3 + ['UNK2']*2)*2 + ['NA+'] + (['UNK1']*3 + ['UNK2']*2) + ['NA+']*4
    assert mpt['resid'] == [1, 1, 1, 2, 2, 3, 3, 3, 4, 4, 5, 6, 6, 6, 7, 7, 8, 9, 10, 11]
    assert mpt['name'] == ['C1', 'C2', 'C3', 'H1', 'N1']*2 + ['NA'] + ['C1', 'C2', 'C3', 'H1', 'N1'] + ['NA']*4
    assert mpt['element'] == ['C', 'C', 'C', 'H', 'N']*2 + ['Na'] + ['C', 'C', 'C', 'H', 'N'] + ['Na']*4
    assert mpt['mass'] == [12, 12, 12, 1, 14]*2 + [23] + [12, 12, 12, 1, 14] + [23]*4
    assert mpt['charge'] == [0.5, -0.5, 0.2, 0, 2]*2 + [1] + [0.5, -0.5, 0.2, 0, 2] + [1]*4

def test_select():
    df1, df2 = getMockTopol()
    
    from mimicpy.io.topol_dict import TopolDict
    topol_dict = TopolDict.from_dict({'MOL1':df1, 'NA1':df2, 'MOL2':df1})
    
    mpt = Mpt([('MOL1', 1), ('NA1', 1), ('MOL2', 1)], topol_dict, 'r')
    
    sele = mpt.select('(type not CT and mol is MOL1) or (type not N* and mol is MOL2)')
    
    assert sele.index.to_list() == [1,2,4,5,7,8,9,10]
    assert sele['type'].to_list() == ['CA']*2+['H', 'N*'] + ['CA']*2 + ['CT', 'H']
    assert sele['resid'].to_list() == [1, 1, 2, 2, 4, 4, 4, 5]
    
    from mimicpy.utils.errors import SelectionError
    
    with pytest.raises(SelectionError) as e:
        assert mpt.select('_type is CA')
    assert str(e.value) == "_type is not a valid selection keyword"
    
    with pytest.raises(SelectionError) as e:
        assert mpt.select('name if CA')
    assert str(e.value) == "if is not a valid logical operator"
    
    with pytest.raises(SelectionError) as e:
        assert mpt.select('name is CA ot name is CT')
    assert str(e.value) == "ot is not a valid boolean operator"
    
    with pytest.raises(SelectionError) as e:
        assert mpt.select('name is CA) or name is CT')
    assert str(e.value) == "Missing open bracket is selection"
    
    with pytest.raises(SelectionError) as e:
        assert mpt.select('name is CA or ( name is CT')
    assert str(e.value) == "Missing closing bracket is selection"
    