""" Module for gro files """

import numpy as np
import pandas as pd
from .parser import Parser


class Gro:
    """ reads gro files """

    def __init__(self, file, mode='r', buffer=1000):
        self.file = file
        self.mode = mode
        self.buffer = buffer
        self._coords = None
        self._box = None
        if mode == 'r':
            self.__read()
        elif mode == 'w':
            pass
        else:  # Raise exception
            pass

    def __read(self):
        """ reads coordinates and box dimensions """

        def mapped(value):
            if value.isnumeric():
                return np.nan
            try:
                return float(value)
            except ValueError:
                return np.nan


        def string_to_array(string):
            return np.array(list(map(mapped, string.split())))

        f = Parser(self.file)
        
        f.readline()
        number_of_atoms = int(f.readline())
        first_atom_line = f.readline()
        f.buffer = len(first_atom_line) * self.buffer

        number_of_rows = len(first_atom_line.split()) - 3
        
        if number_of_rows == 3:
            cols = ['x', 'y', 'z']
        elif number_of_rows == 6:
            cols = ['x', 'y', 'z', 'v_x', 'v_y', 'v_z']
        else:  # Raise exception
            pass

        values = string_to_array(first_atom_line)
        i = 0

        for string in f:
           i += self.buffer
           values = np.append(values, string_to_array(string))

        values = values[~np.isnan(values)]

        expected_len = number_of_atoms * number_of_rows
        if len(values) == expected_len:
            coords = values
            box = f.readline().split()
        elif len(values) == expected_len+3:
            coords = values[:-3]
            box = values[-3:]
        else:  # Raise exception
            pass

        coords = pd.DataFrame(coords.reshape(number_of_atoms, number_of_rows), columns=cols)
        coords['id'] = coords.index.to_numpy()+1

        self._coords = coords.set_index(['id'])
        self._box = box.tolist()

    @property
    def coords(self):
        if self.mode == 'r':
            return self._coords
        self.mode = 'r'
        self.__read()
        return self._coords

    @property
    def box(self):
        if self.mode == 'r':
            return self._box
        self.mode = 'r'
        self.__read()
        return self._box
