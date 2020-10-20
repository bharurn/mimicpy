import time
import pytest
from mimicpy import Gro
from mimicpy.utils.errors import ParserError


def _get_gro_info(gro_file):

    with open(gro_file) as f:
        lines = f.readlines()
        number_of_atoms = int(lines[1].strip())
        box = lines[-1].split()
    
    return number_of_atoms, box

def test_gro1():
    file_to_check = 'gro_files/gro1.gro'

    gro_file = Gro(file_to_check, buffer=1)  # Reading line by line
    coords, box = gro_file.coords, gro_file.box
    number_of_atoms_to_check, box_to_check = _get_gro_info(file_to_check)
    box = [str(a) for a in box]  # Conversion for comparison to naive reader function

    assert box == box_to_check
    assert len(coords) == number_of_atoms_to_check

    gro_file = Gro(file_to_check, buffer=2)  # Reading two lines at a time
    coords, box = gro_file.coords, gro_file.box 
    number_of_atoms_to_check, box_to_check = _get_gro_info(file_to_check)
    box = [str(a) for a in box]  # Conversion for comparison to naive reader function

    assert box == box_to_check
    assert len(coords) == number_of_atoms_to_check

    gro_file = Gro(file_to_check)
    coords, box = gro_file.coords, gro_file.box 
    number_of_atoms_to_check, box_to_check = _get_gro_info(file_to_check)
    box = [str(a) for a in box]  # Conversion for comparison to naive reader function 

    assert set(coords.columns) == {'x', 'y', 'z', 'v_x', 'v_y', 'v_z'}
    assert box == box_to_check
    assert len(coords) == number_of_atoms_to_check

def test_gro2():
    file_to_check = 'gro_files/gro2.gro'

    gro_file = Gro(file_to_check, buffer=1)
    coords, box = gro_file.coords, gro_file.box
    number_of_atoms_to_check, box_to_check = _get_gro_info(file_to_check)
    box = [str(a) for a in box]  # Conversion for comparison to naive reader function

    assert set(coords.columns) == {'x', 'y', 'z'}
    assert box == box_to_check
    assert len(coords) == number_of_atoms_to_check
    
def test_large_gro():
    file_to_check = 'gro_files/4aj3.gro'
    
    start = time.time()
    gro_file = Gro(file_to_check)
    elapsed_time = (time.time() - start) 

    assert elapsed_time < 1.5
    
    coords, box = gro_file.coords, gro_file.box
    number_of_atoms_to_check, box_to_check = _get_gro_info(file_to_check)
    box = [str(a) for a in box]  # Conversion for comparison to naive reader function
    
    assert set(coords.columns) == {'x', 'y', 'z', 'v_x', 'v_y', 'v_z'}
    assert box == box_to_check
    assert len(coords) == number_of_atoms_to_check

def test_bad_gro1():
    file_to_check = 'gro_files/bad_gro1.gro'

    with pytest.raises(ParserError) as error:
        assert Gro(file_to_check)
    assert 'Error parsing gro_files/bad_gro1.gro: Gro file is not formatted properly' in str(error.value)

def test_bad_gro2():
    file_to_check = 'gro_files/bad_gro2.gro'

    with pytest.raises(ParserError) as error:
        assert Gro(file_to_check)
    assert 'Error parsing gro_files/bad_gro2.gro: Gro file is not formatted properly' in str(error.value) 
