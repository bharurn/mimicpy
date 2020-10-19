from mimicpy.core.selector import PyMOL
from mimicpy.utils.errors import MiMiCPyError
import pytest
import pandas as pd
from string import ascii_letters
from random import choice

class MockMPT:
    def __getitem__(self, a):
        df = pd.DataFrame({'name': ['C']*20})
        df.index = df.index+1
        df['id'] = df.index
        return df.loc[a]

def get_mock_xmlrpc_sele():
    dct = {'id':[], 'symbol': [], 'resi_number': [], 'resn': [], 'charge': [], \
            'name': [], 'coord': [], 'type': [], 'mass': [], 'mol': [],\
            'garbage1': [], 'garbage2': [], 'garbage3': []}
    
    with open('prepare_test_files/qmatoms.txt', 'r') as f:
        for i, line in enumerate(f.readlines()):
            if i == 1:
                mol = line.strip()
            elif i > 1:
                name, x, y, z, q = line.split()
                dct['id'].append(i+1)
                dct['type'].append(name)
                dct['name'].append(name)
                dct['symbol'].append(name[0])
                dct['resn'].append(mol)
                dct['mol'].append(mol)
                dct['resi_number'].append(1)
                dct['mass'].append(0)
                dct['coord'].append([float(x), float(y), float(z)])
                dct['charge'].append(float(q))
                dct['garbage1'].append(choice(ascii_letters))
                dct['garbage2'].append(choice(ascii_letters))
                dct['garbage3'].append(choice(ascii_letters))
        
        outer_dct = {'atom': {}}
        outer_dct['atom'] = dct
        return outer_dct
    
class MockPyMOLAtom():
    def __init__(self):
        self.df = pd.DataFrame(get_mock_xmlrpc_sele()['atom'])
        self.i = -1
        
    def __iter__(self):
        return self

    def __next__(self):
        if self.i >= len(self.df)-1:
            raise StopIteration
        else:
            self.i += 1
            return type('obj', (object,), self.df.loc[self.i].to_dict())
    
def test_xmlrpc_pymol(mocker):
    mocker.patch('pymol.cmd.load', return_value=None)
    mocker.patch('pymol.cmd.get_symmetry',\
                 return_value=[40, 40, 40, 90, 90, 90])
    
    mock_sele = get_mock_xmlrpc_sele()
    mocker.patch('pymol.cmd.get_model', return_value=mock_sele)
    
    pym = PyMOL()
     
    assert pym.load(MockMPT(), "--") == [4, 4, 4]
    
    sele = pym.select("--")
    
    assert set(sele.columns) == {"name", "_name", "_element", "_resname", "_resid", 'x', 'y', 'z'}
    
    assert len(sele) == len(mock_sele['atom']['id'])
    
    x,y,z = list(zip(*mock_sele['atom']['coord']))
    assert x == pytest.approx((sele['x']*10).to_list())
    assert y == pytest.approx((sele['y']*10).to_list())
    assert z == pytest.approx((sele['z']*10).to_list())

def test_pymol(mocker):
    mocker.patch('pymol.cmd.load', return_value=None)
    mocker.patch('pymol.cmd.get_symmetry',\
                 return_value=[40, 40, 40, 90, 90, 90])
    
    mocker.patch('pymol.cmd.get_model', return_value=type('obj', (object,), {'atom': MockPyMOLAtom()}))
    
    pym = PyMOL()
     
    assert pym.load(MockMPT(), "--") == [4, 4, 4]
    
    sele = pym.select("--")
    
    assert set(sele.columns) == {"name", "_name", "_element", "_resname", "_resid", 'x', 'y', 'z'}
    
    mock_sele = get_mock_xmlrpc_sele()
    
    assert len(sele) == len(mock_sele['atom']['id'])
    
    x,y,z = list(zip(*mock_sele['atom']['coord']))
    assert x == pytest.approx((sele['x']*10).to_list())
    assert y == pytest.approx((sele['y']*10).to_list())
    assert z == pytest.approx((sele['z']*10).to_list())