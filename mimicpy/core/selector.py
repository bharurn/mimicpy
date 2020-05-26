from ..shell import _local
from .._global import _Global as gbl
from .base import BaseHandle
import xmlrpc.client as xmlrpclib
import pandas as pd

class Selector:
    def loadCoords(self, gro):
        pass
    def select(self, selection):
        pass

class PyMOL(BaseHandle):
    
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
    
    def loadCoords(self, gro):
        if not gbl.host.isLocal() and not self.forceLocal:
            gbl.host.sftp.get(gro, self.downloadpath)
            gro = self.downloadpath
        
        self.cmd.load(gro)
    
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
    
        mpt_sele = self.mpt.selectAtoms(pymol_sele['id'])
        
        # TO DO: check if names/resname, etc. are same and issue warnings accordingly
        
        return mpt_sele.merge(pymol_sele, left_on='id', right_on='id').set_index(['id'])
        