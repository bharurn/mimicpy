from mimicpy.parsers.gro import read
import time

def get_gro_info(gro_file):
    with open(gro_file) as f:
        lines = f.readlines()
        no = int(lines[1].strip())
        box = lines[-1].split()
    
    return no, box

def check_gro(file, gro_data, box, cols):
    gro_data, box = read(file)
    assert set(gro_data.columns) == cols
    
    no, box_ = get_gro_info(file)
    assert all([a == float(b) for a, b in zip(box, box_)] )
    assert len(gro_data) == no

def test_gro_read():
    gro1, box = read('gro_files/gro1.gro')
    check_gro('gro_files/gro1.gro', gro1, box, {'x', 'y', 'z', 'force-x', 'force-y', 'force-z'})
    gro2, box = read('gro_files/gro2.gro')
    check_gro('gro_files/gro2.gro', gro2, box, {'x', 'y', 'z'})
    
def test_gro_time():
    start = time.time()
    gro1, box = read('gro_files/4aj3.gro')
    elapsed_time = (time.time() - start) 
    
    assert elapsed_time < 1.5
    
    check_gro('gro_files/4aj3.gro', gro1, box, {'x', 'y', 'z', 'force-x', 'force-y', 'force-z'})
    
    