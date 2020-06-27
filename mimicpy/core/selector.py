from .._global import _Global as gbl
from ..parsers import gro, parser
from ..utils.errors import MiMiCPyError
import xmlrpc.client as xmlrpclib
import pandas as pd
import tempfile
from collections import defaultdict
from abc import ABC, abstractmethod

###### Selector using Gromacs GRO

class Selector:
    def __init__(self, lines=500):
        self.lines = lines
        
    def load(self, mpt, gro_file):
        self.mpt = mpt
        self.coords, box = gro.read(gro_file, self.lines)
        return box
        
    def select(self, selection):
        """Select MPT atoms and merge with GRO"""
        sele = self.mpt.select(selection)
        return sele.merge(self.coords, left_on='id', right_on='id')

###### Selector using Visualization packages, currently PyMOL and VMD supported
    
class VisPackage(ABC):
    def __init__(self, cmd, load, forcelocal, lines):
        self.cmd = cmd
        self._load = load
        self.forcelocal = forcelocal
        self.lines = lines
    
    def load(self, mpt, gro_file):
        if not gbl.host.isLocal() and self.forceLocal and self.load:
            # if remote and want to run locally, save gro to temp file locally
            temp = tempfile.NamedTemporaryFile(prefix='mimicpy_temp_', suffix=".gro")
            gro_parser = parser.Parser(gro_file, self.lines)
            for i in gro_parser: temp.write(i.encode('utf-8'))
            gro_parser.close()
            gro_file = temp.name
        
        self.mpt = mpt
        if self._load: self._vis_pack_load(gro_file)
        
        return gro.getBox(gro_file)
    
    @abstractmethod
    def _vis_pack_load(self, gro_file):
        pass

    @abstractmethod
    def _sele2df(self, selection):
        """
        Select atoms using the Visualization Package
        and return dataframe, with VisPack columns prefixed with underscore
        Should be implemented in the respecitve VisPack Class
        """
        pass
    
    def select(self, selection):
        sele = self._sele2df(selection)
        mpt_sele = self.mpt[sele['id']]    
        # TO DO: check if names/resname, etc. are same and issue warnings accordingly  
        # the corresp. columns from the vis software will have underscore prefix
        return mpt_sele.merge(sele, left_on='id', right_on='id').set_index(['id'])

class PyMOL(VisPackage):
    """
    PyMOL selector to process a PyMOL selection and combine it with an MPT
    Can be used by connecting to PyMOL using xmlrpc or by executing in the PyMOL interpreter
    
    """
    
    def __init__(self, pymol_handle=None, url='http://localhost:9123', load=True, forcelocal=True, lines=500):
        
        if pymol_handle:
            cmd = pymol_handle
        else:
            try:
                # first try importing pymol in the pymol enviornment
                from pymol import cmd
            except ImportError:
                # if not try connecting by xmlrpc
                cmd = xmlrpclib.ServerProxy(url)
            
                # xmlrpc will silently fail, try a dummy command to check if connection if working
                try:
                    cmd.get_view()
                except ConnectionRefusedError:
                    raise MiMiCPyError(f"Could not connect to PyMOL. Run MiMiCPy in the PyMOL environment or verify that it is running at {url}")
        
        
        super().__init__(cmd, load, forcelocal, lines)
    
    def _vis_pack_load(self, gro_file):
        self.cmd.load(gro_file)
        
    def _sele2df(self, selection):
        sele = self.cmd.get_model(selection, 1)
        
        if isinstance(sele, dict):
            # sele is dict if using xmlrpc
            # atom key of dict has all we need
            df = pd.DataFrame(sele['atom'])
        else:
            # sele will be chempy.models.Indexed object if using from pymol
            # sele.atom is is a list of chempy.Atom objects, each atom has id, symbol, etc.
            df_dict = defaultdict(list)
            params_to_get = ['id', 'symbol', 'name', 'resn', 'resi_number', 'coord']
            for a in sele.atom:
                for i in params_to_get:
                    df_dict[i].append(getattr(a, i))
            df = pd.DataFrame(df_dict, columns = params_to_get)
        
        # extract coordinates and covert from ang to nm
        x,y,z = list(zip(*df[['coord']].apply(lambda x: [i/10 for i in x[0]], axis=1)))
        df.insert(2, "x", x, True) 
        df.insert(2, "y", y, True) 
        df.insert(2, "z", z, True)
        df = df.drop(['coord'], axis=1)
        return df.rename(columns={"name": "_name", "symbol": "_element", "resn": "_resname",\
                               "resi_number": "_resid"})
    
class VMD(VisPackage):
    
    def __init__(self, vmd_handle=None, load=True, forcelocal=True, lines=500):
        if vmd_handle:
            vmd = vmd_handle
        else:
            try:
                import vmd
            except ImportError:
                raise MiMiCPyError("VMD python package not found. Make sure that you have VMD built with python support")
        
        self.molid = 0 # set mol id in case load=False
        
        super().__init__(vmd, load, forcelocal, lines)
    
    def _vis_pack_load(self, gro_file):
        self.molid = self.cmd.molecule.load("gro", gro_file)
    
    def _sele2df(self, selection):
        sele = self.cmd.atomsel(selection, self.molid)
        
        params_to_get = ['name', 'type', 'index', 'mass', 'element', 'resname', 'resid', 'x', 'y', 'z']
        
        df_dict = {}
            
        for i in params_to_get:
            df_dict[i] = getattr(sele, i)
        
        df = pd.DataFrame(df_dict, columns=params_to_get)
        
        df = df.rename(columns={"name": "_name", "element": "_element", "resname": "_resname",\
                               "resid": "_resid", 'mass': '_mass', 'type': '_type'})
        
        return df.rename(columns={"index": "id"})