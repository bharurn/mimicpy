import pandas as pd
from string import ascii_letters
from random import choice
from mimicpy import VMD

class MockMPT:
    def __getitem__(self, a):
        df = get_mock_vmd_sele()
        df.index = df.index+1
        df = df.drop(['garbage1', 'garbage2', 'garbage3', 'x', 'y', 'z', 'index'], axis=1)
        df['id'] = df.index
        return df.loc[a]

def get_mock_vmd_sele():
    dct = {'type': [], 'resid': [], 'resname': [],\
            'name': [], 'charge': [], 'element': [], 'mass': [], 'mol': [], 'x': [], 'y': [], 'z': [],\
            'garbage1': [], 'garbage2': [], 'garbage3': []}
    
    with open('prepare_test_files/qmatoms.txt', 'r') as f:
        for i, line in enumerate(f.readlines()):
            if i == 1:
                mol = line.strip()
            elif i > 1:
                name, x, y, z, q = line.split()
                dct['type'].append(name)
                dct['name'].append(name)
                dct['element'].append(name[0])
                dct['resname'].append(mol)
                dct['mol'].append(mol)
                dct['resid'].append(1)
                dct['mass'].append(0)
                dct['x'].append(float(x))
                dct['y'].append(float(y))
                dct['z'].append(float(z))
                dct['charge'].append(float(q))
                # last three should not be present in the selection (testing clean up of df)
                dct['garbage1'].append(choice(ascii_letters))
                dct['garbage2'].append(choice(ascii_letters))
                dct['garbage3'].append(choice(ascii_letters))
        
        df = pd.DataFrame(dct)
        
        df['index'] = df.index
        return df

def test_vmd(mocker):
    mock_sele = get_mock_vmd_sele()
    
    ## create mock vmd library
    mocker.patch('vmd.molecule.load', return_value=100)
    mocker.patch('vmd.molecule.get_periodic', return_value={'a': 40, 'b': 40, 'c': 40, 'A': 90, 'B': 90, 'C': 90})
    mocker.patch('vmd.atomsel', return_value=mock_sele)
    ###
    
    vmd = VMD(MockMPT(), 'fake.gro') # pass fake gro file so it calls vmd.molecule.load
     
    assert vmd.mm_box == [4, 4, 4]
    
    assert vmd.molid == 100
    
    sele = vmd.select()
    
    assert set(sele.columns) == {"name", "_name", "element", "_element", "resname",\
              "_resname", "resid", "_resid", "mass", '_mass', "charge", "mol", "type", '_type', 'x', 'y', 'z'}
    
    assert len(sele) == len(mock_sele)
    assert (mock_sele['x']/10).to_list() == sele['x'].to_list()
    assert (mock_sele['y']/10).to_list() == sele['y'].to_list()
    assert (mock_sele['z']/10).to_list() == sele['z'].to_list()