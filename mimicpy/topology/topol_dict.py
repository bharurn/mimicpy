"""Module for MiMiCPy-specific molecule:topology dictionary"""


class TopolDict:
    """provides a dictionary with non-repeating topology information"""

    @classmethod
    def from_dict(cls, dict_df):
        keys = list(dict_df.keys())
        repeating = {}
        i = 0
        while i < len(keys):
            key_i = keys[i]
            for j in range(i+1, len(keys)):
                key_j = keys[j]
                if dict_df[key_i].equals(dict_df[key_j]):
                    repeating[key_j] = key_i
                    del dict_df[key_j]
            i += 1
            keys = list(dict_df.keys())
        return cls(dict_df, repeating)

    def __init__(self, dict_df, repeating):
        self.dict_df = dict_df
        self.repeating = repeating

    def __getitem__(self, key):
        if key in self.dict_df:
            return self.dict_df[key]
        if key in self.repeating:
            return self.dict_df[self.repeating[key]]
        raise KeyError(f'Molecule {key} is not in topology.')

    def todict(self):
        extras = self.dict_df.copy()
        for i in self.repeating:
            extras[i] = self.__getitem__(i)
        return extras

    def __str__(self):
        return str(self.todict())

    def __repr__(self):
        return repr(self.todict())

    def keys(self):
        """Handle keys of a TopolDict like keys of a regular dict"""
        combined_dict = self.dict_df.copy()
        combined_dict.update(self.repeating)
        all_keys = combined_dict.keys()
        return all_keys
