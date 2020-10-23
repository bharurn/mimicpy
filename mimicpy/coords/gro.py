"""Module for gro files"""

import numpy as np
import pandas as pd
from ..utils.errors import ParserError
from .base import BaseCoordsClass

class Gro(BaseCoordsClass):
    """reads gro files"""

    def _read(self):
        """Read atom coordinates and box dimensions"""

        def mapped(value):
            if value.isnumeric():
                return np.nan
            try:
                return float(value)
            except ValueError:
                return np.nan

        def string_to_array(string):
            return np.array(list(map(mapped, string.split())))

        self.file.readline()
        number_of_atoms = int(self.file.readline())
        first_atom_line = self.file.readline()
        self.file.buffer = len(first_atom_line) * self.buffer

        number_of_rows = len(first_atom_line.split()) - 3
        if number_of_rows == 3:
            cols = ['x', 'y', 'z']
        elif number_of_rows == 6:
            cols = ['x', 'y', 'z', 'v_x', 'v_y', 'v_z']
        else:
            raise ParserError(self.file_name, details='Gro file is not formatted properly.')

        values = string_to_array(first_atom_line)
        for string in self.file:
            values = np.append(values, string_to_array(string))
        values = values[~np.isnan(values)]

        expected_len = number_of_atoms * number_of_rows
        if len(values) == expected_len:
            coords = values
            box = self.file.readline().split()
        elif len(values) == expected_len + 3:
            coords = values[:-3]
            box = values[-3:]
        else:
            raise ParserError(self.file_name, details='Gro file is not formatted properly.')

        coords = pd.DataFrame(coords.reshape(number_of_atoms, number_of_rows), columns=cols)
        coords['id'] = coords.index.to_numpy()+1

        return coords.set_index(['id']), box.tolist()
        
    def _write(self, mpt_coords):
        pass