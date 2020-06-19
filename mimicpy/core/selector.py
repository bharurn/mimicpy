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
        gro_df, box = gro.read(gro_file, self.lines)
        self.df = mpt.buildSystemTopology()
        self.df = self.df.merge(gro_df, left_on='id', right_on='id')
        return box
        
    def select(self, selection):
        """Translate selection language string into pandas dataframe selection"""
        
        if selection is None:
            #raise SelectionError
            pass
        
        df = self.df # rename df
        
        # if selection is a lambda func, then just call and return it
        # this is provided for debugging puposes
        LAMBDA = lambda:0
        if isinstance(selection, type(LAMBDA)):
            return df[selection(df)]
            
        # below code translates selection langauge to pandas boolean
        # selection eg., resname is SER and id < 25 and mol not Protein_chain_B
        # will be translated to df['resname'] == 'SER' and df.index < 25 and df['mol'] != 'Protein_chain_B'
    
        ev = '' # converted string
        i = 0 # counter to keep track of word position
        for s in selection.split():
            if i == 0: # if starting of set
                ev += f"(df['{s}']"
                # if and/or encountered, reset i to -1 (will become 0 due to i+= 1 at end)
                # so we can start parsing again
            elif s == 'or':
                ev += f' | '
                i = -1
            elif s == 'and':
                ev += f' & '
                i = -1
            elif s == 'is':
                ev += '=='
            elif s == 'not':
                ev += '!='
            elif s == '>' or s == '>=' or s == '<' or s == '<=':
                ev += s
            else: # parse everything else, meant for the third word
                if s.isnumeric():
                    ev += f"{s})"
                else:
                    ev += f"'{s}')"
                
            i += 1

        ev = f"df.loc[{ev}]" # eg., df.loc[ df['resname'] == 'SER' and df.index < 25 and df['mol'] != 'Protein_chain_B' ]
        ev = ev.replace("df['id']","df.index") # replace df['id'] to df.index as id is the index of df
        gbl.logger.write('debug2', f'Selection command translated to: {ev}')
        
        return eval(ev) # evaluate string and return the dataframe

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
    
        mpt_sele = self.mpt.selectByIDs(pymol_sele['id'])
        
        # TO DO: check if names/resname, etc. are same and issue warnings accordingly
        
        return mpt_sele.merge(pymol_sele, left_on='id', right_on='id').set_index(['id'])
        