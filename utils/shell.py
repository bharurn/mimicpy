from . import _local as local, _remote as remote
import os
from stat import S_ISDIR, S_ISREG

class Base():
    
    def __init__(self, directory=".",  modules=[], sources=[], loaders=[], ignloaderr = True):
        
        if directory.strip() == '':
            directory = '.'
        
        ret = self.cd(directory, mkdir=True)
        
        if ret == -1:
            pass
        else:
            if ret == 1:
                print(f"{directory} not found, creating new directory..")
            
            print(f"Setting current directory to {directory}..")
        
        self.modules = modules
        self.sources = sources
        self.loaders = loaders
            
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
        
    def mkdir(self, directory):
        if not self.fileExists(directory):
            os.mkdir(directory)
            return 0
        else:
            return 1
    
    def pwd(self): return os.getcwd()
    
    def read(self, file):
        with open(file, 'r') as f:
                return f.read()
    
    def write(self, content, file):
        with open(file, 'w') as f: f.write(content)
    
    def vi(self, file, mode):
        return open(file, mode)
        
    def run(self, cmd, stdin=None, errorHandle=None, query_rate=0):
        return local.run(cmd, stdin=stdin, errorHandle=errorHandle)
    
    def rename(self, a, b): os.rename(a, b)
    
    def rm(self, a): os.remove(a)
    
    def close(self): pass

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
        
        Base.__init__(self, dir_, modules, sources, loaders)
        
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
    
    def write(self, content, file):
        with self.vi(file, 'w') as f: f.write(content)
        
    def rename(self, a, b): self.sftp.rename(a, b)
    def rm(self, a): self.sftp.remove(a)
    
    def fileExists(self, file):
        try:
            self.sftp.stat(file)
            return True
        except FileNotFoundError:
           return False
    
    def pwd(self): return self.dir
    
    def ls(self, dirc=None, file_eval=lambda a: True, dir_eval=lambda a: True):
        
        files = []
        
        if dirc:
            dirs = self.sftp.listdir_attr(dirc)
        else:
            dirs = self.sftp.listdir_attr()
        
        for entry in dirs:
            if S_ISDIR(entry.st_mode) and dir_eval(entry.filename):
                files.append(entry.filename)
            elif S_ISREG(entry.st_mode) and file_eval(entry.filename):
                files.append(entry.filename)
            
        return files
    
    def vi(self, name, s, prefetch=True):
        file = self.sftp.open(f"{name}", s)
        if prefetch: file.prefetch()
        return file
    
    def sbatch(self, job):
        if job.noCommands():
            raise Exception("No commands in jobscript!")
        job.setDir(self.dir)
        
        f = self.vi(f'{job.name}.sh', 'w')
        f.write(str(job))
        f.close()
        
        def _sbatch_err(err):
            if 'error' in err.lower():
                raise Exception(err) 
            
        out = self.run(f'cd {self.dir} && sbatch {job.name}.sh', errorHandle=_sbatch_err, onNewChan=True)
        
        try:
            idx = int(out.split()[3])
        except:
            raise Exception(out)
        
        return idx
        
    def cd(self, directory, mkdir=False):
        if directory.strip() == '.':
            return -1
        
        if not self.fileExists(directory):
             if mkdir:       
                self.sftp.mkdir(directory)
                return 1
             else:
                 raise Exception(f'Directory {directory} not found')
        else:     
            self.dir = directory+'/'
            self.run(f'cd {self.dir}')
            self.sftp.chdir(self.dir)
            return 0
    
    def mkdir(self, directory):
        if not self.fileExists(directory):
            self.sftp.mkdir(directory)
            return 0
        else:
            return 1
    
    def cp(self, f1, f2):
        self.run(f"cp {f1} {f2}")
    
    def scancel(self, jobid):
        return self.run(f'scancel {jobid}', onNewChan=True)
    
    def close(self): self.__del__()