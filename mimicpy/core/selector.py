from .._global import _Global as gbl
from ..parsers import gro, parser
import xmlrpc.client as xmlrpclib
import pandas as pd
import tempfile
from collections import defaultdict

class Selector:
    def __init__(self, lines=500):
        self.lines = lines
        
    def load(self, mpt, gro_file):
        self.mpt = mpt
        self.gro_df, box = gro.read(gro_file, self.lines)
        return box
        
    def select(self, selection):
        """Select MPT atoms and merge with GRO"""
        atoms = self.mpt.select(selection)
        atoms = atoms.merge(self.gro_df, left_on='id', right_on='id')
        return atoms

class PyMOL(Selector):
    """PyMOL selector to process a PyMOL selection and combine it with an MPT
    Can be used by connecting to PyMOL using xmlrpc or by executing in the PyMOL interpreter
    """
    
    def __init__(self, cmd=None, url='http://localhost:9123', load=True, forcelocal=False, lines=500):
        
        if cmd is None:
            self.cmd = xmlrpclib.ServerProxy(url)
        else:
            self.cmd = cmd
        
        self.load = load
        self.forcelocal = forcelocal
        self.lines = lines
    
    def load(self, mpt, gro):
        if not gbl.host.isLocal() and self.forceLocal:
            # if remote and want to run pymol locally, download gro to temp file
            temp = tempfile.NamedTemporaryFile(prefix='mimicpy_temp_', suffix=".gro")
            gro_parser = parser.Parser(gro, self.lines)
            for i in gro_parser: temp.write(i.encode('utf-8'))
            gro_parser.close()
            gro = temp.name
        
        self.cmd.load(gro)
        self.mpt = mpt
        
        return gro.getBox(gro)
    
    def __sele2df(self, sele):
        if isinstance(sele, dict):
            # sele is dict if using xmlrpc
            # atom key of dict has all we need
            df_dict = sele['atom']
        else:
            # sele will br chempy.models.Indexed object if using from pymol
            # sele.atom is is a list of chempy.Atom objects, each atom has id, symbol, etc.
            df_dict = defaultdict(list)
            params_to_get = ['id', 'symbol', 'name', 'resn', 'resi_number', 'coord']
            for a in sele.atom:
                for i in params_to_get:
                    df_dict[i] = getattr(a, i)
        
        df = pd.DataFrame(df_dict)
        # extract coordinates and covert from ang to nm
        x,y,z = list(zip(*df[['coord']].apply(lambda x: [i/10 for i in x[0]], axis=1)))
        df.insert(2, "x", x, True) 
        df.insert(2, "y", y, True) 
        df.insert(2, "z", z, True)
        df = df.drop(['coord'], axis=1)
        df = df.rename(columns={"name": "pm_name", "symbol": "pm_symbol", "resn": "pm_resn",\
                               "resi_number": "pm_resi_number"})
    
    def select(self, selection=None):
        if selection is None: selection = 'sele'
        sele_return = self.cmd.get_model(selection, 1)
        pymol_sele = self.__sele2df(sele_return)
    
        mpt_sele = self.mpt[pymol_sele['id']]
        
        # TO DO: check if names/resname, etc. are same and issue warnings accordingly
        
        return mpt_sele.merge(pymol_sele, left_on='id', right_on='id').set_index(['id'])
        