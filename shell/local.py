from pyshell import local
from pygmx.shell import base
import os

class Local(base.Base):
    
    def __init__(self, modules=[], sources=[], directory=".", loaders=[], ignloaderr = True):
        super().__init__(modules, sources, directory, loaders, ignloaderr)
        
        self.name = 'localhost'
    
    def cd(self, directory, mkdir=False):
        return local.cd(directory, mkdir)
    
    def checkFile(self, file, throw=False):
        ret = local.fileExists(file)
        
        if ret == True:
            print(f"{file} found!")
        elif throw:
            raise Exception(f"{file} not found!")
    
    def fileExists(self, file):
        local.fileExists(file)
    
    def pwd(): return os.getcwd()
    
    def read(self, file):
        with open(file, 'r') as f:
                return f.read()
    
    def write(self, content, file, mode='w'):
        with open(file, mode) as f: f.write(content)
        
    def run(self, cmd, stdin=None, errorHandle=None):
        return local.run(cmd, stdin=stdin, errorHandle=errorHandle)
    
    def rename(self, a, b): os.rename(a, b)
    
    def rm(self, a): os.remove(a)