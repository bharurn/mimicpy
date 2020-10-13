import logging
import pandas as pd
from ..io.mpt import Mpt
from ..io.gro import Gro
from ..utils.errors import MiMiCPyError
from ._tclvmd import TclVMDConnector
import xmlrpc.client as xmlrpclib
from collections import defaultdict
from abc import ABC, abstractmethod


class GroSelector:

    def __init__(self, mpt_file, gro_file, buffer=1000):
        self.mpt = Mpt.from_file(mpt_file, buffer=buffer)
        self.gro = Gro(gro_file, buffer=buffer)
        if self.mpt.number_of_atoms != len(self.gro.coords):
            raise MiMiCPyError("Number of atoms in mpt and number of atoms in gro do not match.")

    @property
    def mm_box(self):
        return self.gro.box

    def select(self, selection):
        """Select MPT atoms and merge with GRO"""
        sele = self.mpt.select(selection)
        df = sele.merge(self.gro.coords, left_on='id', right_on='id')

        if df.empty:
            raise MiMiCPyError("The atoms selected from mpt were not found in gro file")
        return df

###### Selector using Visualization packages, currently PyMOL and VMD supported

class VisPackage(ABC):
    ######Core Methods
    ##
    def __init__(self, mpt_file, gro_file, cmd):
        self.cmd = cmd
        self.mpt = Mpt.from_file(mpt_file)
        if gro_file:
            self._vis_pack_load(gro_file)

    def select(self, selection):
        sele = self._sele2df(selection)
        mpt_sele = self.mpt[sele['id']]
        # TO DO: check if names/resname, etc. are same and issue warnings accordingly
        # the corresp. columns from the vis software will have underscore prefix
        df = mpt_sele.merge(sele, left_on='id', right_on='id').set_index(['id'])

        if df.empty:
            raise MiMiCPyError("The atoms IDs in selected from the visualization software do not exist in the mpt file")
        return df
    ##
    ######

    ######Methods to override in children
    ##
    @abstractmethod
    def _vis_pack_load(self, gro_file):
        # call load gro file of vis pack
        pass

    @abstractmethod
    def get_box(self):
        # return box size is in nm
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

    def __init__(self, mpt_file, gro_file, url=None):

        if url is None:
            try:
                # first try importing pymol in the pymol enviornment
                from pymol import cmd
            except ImportError:
                raise MiMiCPyError(f"Could not connect to PyMOL, make sure PyMOL is installed")
        else:
            # if not try connecting by xmlrpc
            cmd = xmlrpclib.ServerProxy(url)

            # xmlrpc will silently fail, try a dummy command to check if connection if working
            try:
                cmd.get_view()
            except ConnectionRefusedError:
                raise MiMiCPyError(f"Could not connect to PyMOL xmlrpc server at address {url}")


        super().__init__(mpt_file, gro_file, cmd)

    def _vis_pack_load(self, gro_file):
        self.cmd.load(gro_file)

    def get_box(self):
        box = self.cmd.get_symmetry("all")
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

    def __init__(self, mpt_file, gro_file, tcl_vmd_params=None):
        if tcl_vmd_params:
            # use the Tcl connector
            vmd = TclVMDConnector(tcl_vmd_params)
        else:
            try:
                import vmd
            except ImportError:
                raise MiMiCPyError("VMD python package not found. Make sure that you have VMD built with python support"
                                   "\nConsider using the Tcl VMD Connector instead.")

        self.molid = -1 # default to top mol

        super().__init__(mpt_file, gro_file, vmd)

    def _vis_pack_load(self, gro_file):
        self.molid = self.cmd.molecule.load("gro", gro_file)

    def get_box(self):
        box = self.cmd.molecule.get_periodic(self.molid)
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
                df_dict[i] = [j+1 for j in df_dict[i]]

            if i in ['x', 'y', 'z']:
                df_dict[i] = [j/10 for j in df_dict[i]]

        df = pd.DataFrame(df_dict, columns=params_to_get)

        df = df.rename(columns={"name": "_name", "element": "_element", "resname": "_resname",\
                               "resid": "_resid", 'mass': '_mass', 'type': '_type'})

        return df.rename(columns={"index": "id"})
