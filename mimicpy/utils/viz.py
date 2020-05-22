from ..shell import _local
from .._global import _Global as gbl
import xmlrpc.client as xmlrpclib

class PyMol:
    
    def _hook(self, out):
        if 'xml-rpc server running on host' in out:
            
            self.port = out.split()[-1]
            self.host = out.split()[-3].replace(',', '')
            
            return True
        else:
            return False
    
    def connect(self, launch=True, host='localhost', port=9123):
        if launch:
            _local.runbg('pymol -R', hook=self._hook)
        else:
            self.host = host
            self.port = port
        
        try:
            self.cmd = xmlrpclib.ServerProxy(f'http://{host}:{port}')
        except:
            print('Cannot connect to PyMol at host {host} and port {port}!')
    
    def loadCoords(self, gro, forceLocal=False, downloadTo='temp.gro'):
        if not gbl.host.isLocal and not forceLocal:
            gbl.host.sftp.get(gro, downloadTo)
            gro = downloadTo
        
        self.cmd.load(gro)
        