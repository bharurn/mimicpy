from pyshell import remote
from pygmx.shell import base

class SSH(remote.shell.Shell, base.Base):
    
    def __init__(self, server, modules=[], sources=[], directory=".", loaders=[], ignloaderr = True):
        
        remote.shell.Shell.__init__(self, server)
        
        base.Base.__init__(self, modules, sources, directory, loaders, ignloaderr)
        
        self.query_rate = 3
    
    def checkFile(self, file, throw=False):
        ret = self.fileExists(file)
        
        if ret == True:
            print(f"{file} found!")
        elif throw:
            raise Exception(f"{file} not found!")
        else:
            return False
    
    def read(self, file):
        with self.vi(file, 'r') as f:
            return f.read().decode('utf-8')
    
    def write(self, content, file, mode='w'):
        with self.vi(file, mode) as f: f.write(content)
        
    def rename(self, a, b): self.sftp.rename(a, b)
    def rm(self, a): self.sftp.remove(a)