import mimicpy
import pickle

def test_dppc():
    log = mimicpy.utils.logger.LogString()
    mimicpy.setLogger(3, log)
    
    warns = mimicpy.utils.logger.LogString()
    mimicpy.redirectWarnings(warns)
    
    mol_list, topol_dict = mimicpy.parsers.top.read('dppc/topol.top', guess_elems=False)
    
    files = ['topol.top', 'gromos53a6.ff/forcefield.itp', 'gromos53a6.ff/ffnonbonded.itp', 'dppc.itp', 'gromos53a6.ff/spc.itp']
    files = ['dppc/'+i for i in files]
    
    for i in log.splitlines():
        assert i.split()[0] in files
    
    assert warns == 'Cannot find dppc/posre.itp, skipping..\n'
    
    assert mol_list==[('DPPC', 2), ('SOL', 100)]
    
    with open('dppc/result.pkl', 'rb') as f:
        known_topol_dict = pickle.load(f)
    
    assert topol_dict.repeating == known_topol_dict.repeating
    
    assert topol_dict.dict_df.keys() == known_topol_dict.dict_df.keys()
    
    for k,v in topol_dict.dict_df.items():
        assert v.equals(v)

def test_4aj3():
    log = mimicpy.utils.logger.LogString()
    mimicpy.setLogger(3, log)
    
    warns = mimicpy.utils.logger.LogString()
    mimicpy.redirectWarnings(warns)
    
    mol_list, topol_dict = mimicpy.parsers.top.read('4aj3/topol.top', guess_elems=True)
    
    assert mol_list==[('Protein', 1), ('NAP', 1), ('ICT', 1), ('SOL', 47708), ('NA', 18), ('SOL', 18)]