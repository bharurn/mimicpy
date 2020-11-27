from abc import ABC, abstractmethod
from collections import defaultdict
import logging
import xmlrpc.client as xmlrpclib
import pandas as pd
from ..topology.mpt import Mpt
from ..coords.base import CoordsIO
from ..utils.errors import MiMiCPyError
from ..utils.strings import print_table


class DefaultSelector:

    def __init__(self, mpt_file, coord_file, buffer=1000, nonstandard_atomtypes=None, gmxdata=None, file_ext=None):
        self.mpt = Mpt.from_file(mpt_file, buffer=buffer, nonstandard_atomtypes=nonstandard_atomtypes,\
                                 gmxdata=gmxdata, file_ext=file_ext)
        self.coords_reader = CoordsIO(coord_file, buffer=buffer)
        n_mpt = self.mpt.number_of_atoms
        n_coords = len(self.coords_reader.coords)
        if n_mpt != n_coords:
            raise MiMiCPyError('Number of atoms in topology and coordinates do not match ({} vs {})'.format(n_mpt, n_coords))

    @property
    def mm_box(self):
        return self.coords_reader.box

    def select(self, selection):
        """Select MPT atoms and merge with GRO"""
        sele = self.mpt.select(selection)
        df = sele.merge(self.coords_reader.coords, left_on='id', right_on='id')

        if df.empty:
            raise MiMiCPyError('The atoms selected from mpt were not found in the coordinate file')
        return df

###### Selector using Visualization packages, currently PyMOL and VMD supported

class VisPackage(ABC, DefaultSelector):
    ######Core Methods
    ##
    def __init__(self, mpt_file, coord_file, cmd, buffer, nonstandard_atomtypes, gmxdata, file_ext):
        self.cmd = cmd
        self.mpt = Mpt.from_file(mpt_file, buffer=buffer, nonstandard_atomtypes=nonstandard_atomtypes,\
                                 gmxdata=gmxdata, file_ext=file_ext)
        if coord_file:
            self._vis_pack_load(coord_file)

    def select(self, selection=None):
        sele = self._sele2df(selection)
        mpt_sele = self.mpt[sele['id']]
        df = mpt_sele.merge(sele, left_on='id', right_on='id').set_index(['id'])
        
        if df.empty:
            raise MiMiCPyError('The atoms IDs in selected do not exist in {}'.format(self.mpt))
            
        # check for mismatch between mpt and mol viz
        cols = [i for i in df.columns if i.startswith('_')]
        for j in cols:
            i = j[1:]
            logging.warning("\nThe following atom(s) did that have matching '{}' information:\n".format(i))
            d = df[df[i] != df[j]]
            dct = {'Atom ID': d.index.to_list(), 'From Topology':d[i].to_list(), 'From Software': d[j].to_list()}
            print_table(dct, logging.warning)

        return df

    ##
    ######

    ######Methods to override in children
    ##
    @abstractmethod
    def _vis_pack_load(self, coord_file):
        """
        Call the corrseponding coord load function of Visualization Package
        """
        pass

    @property
    @abstractmethod
    def mm_box(self):
        """
        Call the corrseponding function of Visualization Package to get box size
        """
        pass

    @abstractmethod
    def _sele2df(self, selection):
        """
        Select atoms using the Visualization Package
        and return dataframe, with VisPack columns prefixed with underscore
        Should be implemented in the respecitve VisPack Class
        """
        pass
    ##
    #######

class PyMOL(VisPackage):
    """
    PyMOL selector to process a PyMOL selection and combine it with an MPT
    Can be used by connecting to PyMOL using xmlrpc or by executing in the PyMOL interpreter
    """

    def __init__(self, mpt_file, coord_file=None, url=None, buffer=1000, nonstandard_atomtypes=None, gmxdata=None, file_ext=None):

        if url is None:
            try:
                # import pymol in the pymol enviornment
                from pymol import cmd
            except ImportError:
                raise MiMiCPyError('Could not connect to PyMOL, make sure PyMOL is installed')
        else:
            # connecting by xmlrpc url
            cmd = xmlrpclib.ServerProxy(url)

            # xmlrpc will silently fail, try a dummy command to check if connection is working
            try:
                cmd.get_view()
            except ConnectionRefusedError:
                raise MiMiCPyError('Could not connect to PyMOL xmlrpc server at address {}'.format(url))


        super().__init__(mpt_file, coord_file, cmd, buffer, nonstandard_atomtypes, gmxdata, file_ext)

    def _vis_pack_load(self, coord_file):
        self.cmd.load(coord_file)

    @property
    def mm_box(self):
        box = self.cmd.get_symmetry("all")
        # convert from ang to nm
        return [b/10 for b in box[:3]]


    def _sele2df(self, selection):
        if selection is None:
            selection = 'sele'

        sele = self.cmd.get_model(selection, 1)

        params_to_get = ['id', 'symbol', 'name', 'resn', 'resi_number', 'coord']

        if isinstance(sele, dict):
            # sele is dict if using xmlrpc
            # atom key of dict has all we need
            df_dict = sele['atom']
        else:
            # sele will be chempy.models.Indexed object if using from pymol
            # sele.atom is is a list of chempy.Atom objects, each atom has id, symbol, etc.
            df_dict = defaultdict(list)
            for a in sele.atom:
                for i in params_to_get:
                    df_dict[i].append(getattr(a, i))

        df = pd.DataFrame(df_dict, columns=params_to_get)

        # extract coordinates and covert from ang to nm
        x, y, z = list(zip(*df[['coord']].apply(lambda x: [i/10 for i in x[0]], axis=1)))
        df.insert(2, "x", x, True)
        df.insert(2, "y", y, True)
        df.insert(2, "z", z, True)
        df = df.drop(['coord'], axis=1)
        return df.rename(columns={"name": "_name", "symbol": "_element", "resn": "_resname",\
                               "resi_number": "_resid"})

class VMD(VisPackage):

    def __init__(self, mpt_file, coord_file=None, buffer=1000, nonstandard_atomtypes=None, gmxdata=None, file_ext=None):
        try:
            import vmd
        except ImportError:
            raise MiMiCPyError("VMD python package not found. Make sure that you have VMD built with python support"
                                   "\nConsider using the Tcl VMD Connector instead.")

        self.molid = -1 # default to top mol

        super().__init__(mpt_file, coord_file, vmd, buffer, nonstandard_atomtypes, gmxdata, file_ext)

    def _vis_pack_load(self, coord_file):
        self.molid = self.cmd.molecule.load(coord_file.split('.')[-1], coord_file)

    @property
    def mm_box(self):
        box = self.cmd.molecule.get_periodic(self.molid)
        # convert from ang to nm
        return [box[k]/10 for k in ['a', 'b', 'c']]

    def _sele2df(self, selection):
        if selection is None:
            selection = 'atomselect0'

        sele = self.cmd.atomsel(selection, self.molid)

        params_to_get = ['name', 'type', 'index', 'mass', 'element', 'resname', 'resid', 'x', 'y', 'z']

        df_dict = {}

        for i in params_to_get:
            df_dict[i] = getattr(sele, i)

            if i == 'index':
                # vmd uses 0 based index, mimicpy uses 1 based index
                df_dict[i] = [j+1 for j in df_dict[i]]

            if i in ['x', 'y', 'z']:
                # convert from ang to nm
                df_dict[i] = [j/10 for j in df_dict[i]]

        df = pd.DataFrame(df_dict, columns=params_to_get)

        df = df.rename(columns={"name": "_name", "element": "_element", "resname": "_resname",\
                               "resid": "_resid", 'mass': '_mass', 'type': '_type'})

        return df.rename(columns={"index": "id"})
