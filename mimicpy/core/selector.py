from ..shell import _local
from .._global import _Global as gbl
from ..parsers import gro
import xmlrpc.client as xmlrpclib
import pandas as pd

class Selector:
    def __init__(self, lines=500):
        self.lines = lines
        
    def load(self, mpt, gro_file):
        gro_df, box = gro.read(gro_file, self.lines)
        self.df = mpt.getDF()
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
    
    def _hook(self, out):
        if 'xml-rpc server running on host' in out:
            
            self.port = out.split()[-1]
            self.host = out.split()[-3].replace(',', '')
            
            return True
        else:
            return False
    
    def __init__(self, launch=True, host='localhost', port=9123, load=True, forcelocal=False, downloadpath='temp.gro'):
        if launch:
            _local.runbg('pymol -R', hook=self._hook)
        else:
            self.host = host
            self.port = port
        
        try:
            self.cmd = xmlrpclib.ServerProxy(f'http://{host}:{port}')
        except:
            print('Cannot connect to PyMol at host {host} and port {port}!')
        
        self.load = load
        self.forcelocal = forcelocal
        self.downloadpath = downloadpath
    
    def load(self, mpt, gro):
        if not gbl.host.isLocal() and not self.forceLocal:
            gbl.host.sftp.get(gro, self.downloadpath)
            gro = self.downloadpath
        
        self.cmd.load(gro)
        self.mpt = mpt
        
        return gro.getBox(gro)
    
    def select(self, selection=None):
        if selection is None: selection = 'sele'
        ids = self.cmd.get_model(selection, 1)
        pymol_sele = pd.DataFrame(ids['atom'])
        # extract coordinates and covert from ang to nm
        x,y,z = list(zip(*pymol_sele[['coord']].apply(lambda x: [i/10 for i in x[0]], axis=1)))
        pymol_sele.insert(2, "x", x, True) 
        pymol_sele.insert(2, "y", y, True) 
        pymol_sele.insert(2, "z", z, True)
        pymol_sele = pymol_sele.drop(['coord'], axis=1)
        pymol_sele = pymol_sele.rename(columns={"name": "pm_name", "symbol": "pm_symbol", "resn": "pm_resn",\
                               "resi_number": "pm_resi_number"})
    
        mpt_sele = self.mpt.selectByIDs(pymol_sele['id'])
        
        # TO DO: check if names/resname, etc. are same and issue warnings accordingly
        
        return mpt_sele.merge(pymol_sele, left_on='id', right_on='id').set_index(['id'])
        