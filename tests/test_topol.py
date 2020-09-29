import mimicpy
from mimicpy.utils.errors import ParserError
from mimicpy.io.top import Top
import pickle
import pytest

def test_dppc():
#    log = mimicpy.utils.logger.LogString()
#    mimicpy.setLogger(3, log)

#    warns = mimicpy.utils.logger.LogString()
#    mimicpy.redirectWarnings(warns)

    top = Top('dppc/topol.top')
#    mol_list, topol_dict = top.read()
    mol_list = top.molecules
    topol_dict = top.topol_dict
#    mol_list, topol_dict = mimicpy.parsers.top.read('dppc/topol.top', guess_elems=False)

    files = ['topol.top', 'gromos53a6.ff/forcefield.itp', 'gromos53a6.ff/ffnonbonded.itp', 'dppc.itp', 'gromos53a6.ff/spc.itp']
    files = ['dppc/'+i for i in files]

#    for i in log.splitlines():
#        assert i.split()[0] in files

 #   assert warns == 'Cannot find dppc/posre.itp, skipping..\n'

    assert mol_list==[('DPPC', 2), ('SOL', 100)]

#    with open('dppc/result.pkl', 'rb') as f:
#        known_topol_dict = pickle.load(f) # Still uses the old parsers module??

#    assert topol_dict.repeating == known_topol_dict.repeating

#    assert topol_dict.dict_df.keys() == known_topol_dict.dict_df.keys()

#    for k,v in topol_dict.dict_df.items():
#        assert v.equals(v)

def test_4aj3():
#    log = mimicpy.utils.logger.LogString()
#    mimicpy.setLogger(3, log)

#    warns = mimicpy.utils.logger.LogString()
#    mimicpy.redirectWarnings(warns)

#    with pytest.raises(ParserError) as e:
#        assert mimicpy.io.top.read('4aj3/topol.top', guess_elems=False)
#    assert "Cannot determine atomic symbol for atom with name C6N and type NAP_CA in residue NAP" in str(e.value)

    top = Top('4aj3/topol.top', mode='w')
    mol_list, topol_dict = top.molecules, top.topol_dict

    assert mol_list==[('Protein', 1), ('NAP', 1), ('ICT', 1), ('SOL', 47708), ('NA', 18), ('SOL', 18)]
    assert set(topol_dict.keys()) == {'Protein', 'NAP', 'ICT', 'SOL', 'NA', 'SOL'}
    assert set(topol_dict['NAP']['element']) == {'C', 'H', 'O', 'N', 'P'}
    assert set(topol_dict['ICT']['element']) == {'C', 'H', 'O'}