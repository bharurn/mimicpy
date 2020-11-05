from mimicpy.utils.errors import ParserError
from mimicpy import Top
import pickle
import pytest
import io, logging

def test_dppc():
    warns = io.StringIO()

    fileh = logging.StreamHandler(warns)
    formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
    fileh.setFormatter(formatter)

    log = logging.getLogger()  # root logger
    for hdlr in log.handlers[:]:  # remove all old handlers
        log.removeHandler(hdlr)
    log.addHandler(fileh)      # set the new handler

    top = Top('dppc/topol.top')
    mol_list = top.molecules
    topol_dict = top.topol_dict

    files1 = ['topol.top', 'forcefield.itp']
    files2 = ['dppc_A.itp', 'dppc_B.itp', 'spc.itp']

    for i in files1:
        assert 'No atoms found in {}'.format(i) in warns.getvalue()

    for i in files2:
        assert 'Read atoms from {}'.format(i) in warns.getvalue()

    assert 'Could not find posre.itp in local or Gromacs data directory. Skipping...' in warns.getvalue()

    assert mol_list==[('DPPC_A', 1), ('DPPC_B', 2), ('SOL', 100)]

    with open('dppc/result.pkl', 'rb') as f:
        known_topol_dict = pickle.load(f) # Still uses the old parsers module??

    assert topol_dict.repeating == {'DPPC_B': 'DPPC_A'}

    assert topol_dict.dict_df.keys() == known_topol_dict.keys()

    for k,v in topol_dict.dict_df.items():
        assert v.equals(v)

def test_4aj3():
    with pytest.raises(ParserError) as e:
        assert Top('4aj3/topol.top', guess_elements=False)
    assert "Cannot determine atomic symbol for atom with name C6N and type NAP_CA in residue NAP" in str(e.value)

    top = Top('4aj3/topol.top', mode='w')
    mol_list, topol_dict = top.molecules, top.topol_dict

    assert mol_list==[('Protein', 1), ('NAP', 1), ('ICT', 1), ('SOL', 47708), ('NA', 18), ('SOL', 18)]
    assert set(topol_dict.keys()) == {'Protein', 'NAP', 'ICT', 'SOL', 'NA', 'SOL'}
    assert set(topol_dict['NAP']['element']) == {'C', 'H', 'O', 'N', 'P'}
    assert set(topol_dict['ICT']['element']) == {'C', 'H', 'O'}