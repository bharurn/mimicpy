from . import _local as local, _remote as remote
import os

class Base():
    
    def __init__(self, directory="."):
        
        if directory.strip() == '':
            directory = '.'
        
        ret = self.cd(directory, mkdir=True)
        
        if ret == -1:
            pass
        else:
            if ret == 1:
                print(f"{directory} not found, creating new directory..")
            
            print(f"Setting current directory to {directory}..")
            
class Local(Base):
    
    def __init__(self, directory="."):
        super().__init__(directory)
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

class Remote(remote.Shell, Base):
    
    def __init__(self, work_dir, ssh_config, modules=[], sources=[], loaders=[], ignloaderr = True):
        
        if ':' in work_dir:
            r = work_dir.split(':')
            server = r[0]
            dir_ = r[1]
            if dir_ == '': dir_ = '.'
        else:
            server = work_dir
            dir_ = '.'
            
        print(f"Setting remote machine {r[0]} as host..")
        
        remote.Shell.__init__(self, server, ssh_config)
        
        Base.__init__(self, dir_)
        
        self.query_rate = 3
         
        for module in modules:
            print(f"Loading module {module}..")
            try:
                self.run(f'module load {module}')
            except Exception as e:
                if ignloaderr:
                    print(f"Warning! {e}. Ignoring..")
                else:
                    raise Exception(str(e))
        
        for source in sources:
            print(f"Sourcing {source}..")
            try:
                self.run(f'source {source}')
            except Exception as e:
                if ignloaderr:
                    print(f"Warning! {e}. Ignoring..")
                else:
                    raise Exception(str(e))
        
        for loader in loaders:
            print(f"Running commands: {loader}..")
            try:
                self.run(loader)
            except Exception as e:
                if ignloaderr:
                    print(f"Warning! {e}. Ignoring..")
                else:
                    raise Exception(str(e))
      
    def source(self, sources):
         for source in sources:
            self.run(f'source {source}')
        
    def module(self, modules):
        for module in modules:
            self.run(f'module load {module}')
    
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
    
    def close(self): self.__del__()