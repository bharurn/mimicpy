from . import _local as local, _remote as remote
import os
from stat import S_ISDIR, S_ISREG
from shutil import copyfile
import re
from ..utils.errors import SlurmBatchError
from ..utils.fstring import f
from .._global import _Global as _gbl

class Base():
    
    def __init__(self, directory,  *loaders):
        
        if directory.strip() == '':
            directory = '.'
        
        ret = self.cd(directory, mkdir=True)
        
        if ret == -1:
            pass
        else:
            if ret == 1:
                _gbl.logger.write('debug', f"{directory} not found, creating new directory..")
            
            _gbl.logger.write('debug', f"Setting current directory to {directory}..")
        
        self.loaders = []
        self.loader_str = ''
        self.loader_out = ''
        
        if loaders:
            self.addLoaders(*loaders)
        
    def addLoaders(self, *loader):
         self.loaders.extend(loader)
         if self.loader_str.strip() == '':
             self.loader_str += ' ; '.join(loader)
         else:
             self.loader_str = ' ; '.join(loader)
        
         print_hook = lambda cmd,out: _gbl.logger.write('debug2', f'Running loader {cmd}..\n{out}')
        
         self.loader_out = self.run(self.loader_str, hook = print_hook, fresh=True)
         
    def rename(self, a, b): self.hndl().rename(a, b)
    def rm(self, a): self.hndl().remove(a)
    def pwd(self): return self.hndl().getcwd()
    
    def checkFile(self, file, throw=False):
        ret = self.fileExists(file)
        
        if ret == True:
            _gbl.logger.write('debug', f"{file} found!")
        elif throw:
            raise FileNotFoundError(f"{file} not found!")
    
    def read(self, file):
        with self.open(file, 'r') as f:
            out = f.read()
            try:
                return out.decode()
            except (UnicodeDecodeError, AttributeError):
                return out
    
    def write(self, content, file):
        with self.open(file, 'w') as f: f.write(content)
        
    def mkdir(self, directory):
        if not self.fileExists(directory):
            self.hndl().mkdir(directory)
            return 0
        else: return 1
        
    def cd(self, directory, mkdir=False):
        if directory.strip() == '.':
            return -1
        
        if not self.fileExists(directory):
             if mkdir:       
                self.hndl().mkdir(directory)
                return 1
             else:
                 raise FileNotFoundError(f'Directory {directory} not found')
        else:
            self.hndl().chdir(directory+'/')
            return 0
    
    def run(self, cmd, stdin=None, hook=None, fresh=False, dirc='', combine=False):
        
        if not fresh:
            cmd = self.loader_str + ' ; ' + cmd
        
        out, err = self._run(cmd, stdin=stdin, dirc=dirc)
           
        if not fresh:
            out = out.replace(self.loader_out, '')
            out = out.replace(self.loader_err, '')
            
            err = err.replace(self.loader_err, '')
            err = err.replace(self.loader_out, '')
    
        if hook: hook(err+'\n'+out)
        
        if combine:
            return err+'\n'+out
        else:
            return out, err
    
    def sbatch(self, job, dirc=''):
        jbs = f'{dirc}/{job.name}.sh'
        
        self.write(str(job), jbs)
        
        if job.noCommands():
            raise SlurmBatchError(jbs, "No commands found!")
        
        jid = 0
        def _sbatch_err(txt):
            if 'error' in txt.lower():
                raise SlurmBatchError(jbs, txt)
            else:
                match = re.search(r'Submitted batch job (\w*)', txt)
                if match:
                    global jid
                    jid = match.groups()[0]
                else: raise SlurmBatchError(jbs, txt)
        
        self.run(f'sbatch {job.name}.sh', hook=_sbatch_err, dirc=dirc, fresh=True)
        
        return jid
    
    def scancel(self, jobid):
        return self.run(f'scancel {jobid}', fresh=True)
        
class Local(Base):
    
    def __init__(self, directory, shell_path, *loaders):
        self.name = 'localhost'
        if shell_path:
            self.shell_path = shell_path
        else:
            self.shell_path = os.environ['SHELL'] # check env vars and get shell exec, works only for UNIX
        super().__init__(directory, *loaders)
    
    def hndl(self): return os
    
    def fileExists(self, file):
        return os.path.isfile(file) or os.path.isdir(file)
    
    def open(self, file, mode):
        if 'r' in mode: self.checkFile(file, throw=True)
        return open(file, mode)
        
    def run(self, cmd, stdin=None, hook=None, fresh=False, dirc=''):
        replace = ''
        if not fresh:
            cmd = self.loader_str + ' ; ' + cmd
            replace = self.loader_out
        if dirc != '': cmd = f'cd {dirc} ; ' + cmd
        return local.run(cmd, self.shell_path, replace, stdin=stdin, hook=hook)
    
    def runbg(self, cmd, hook=None, fresh=False, dirc='', query_rate=1):
        replace = ''
        if not fresh:
            cmd = self.loader_str + ' ; ' + cmd
            replace = self.loader_out
        if not dirc != '': cmd = f'cd {dirc} ; ' + cmd
        return local.runbg(cmd, self.shell_path, replace, hook=hook)

    def ls(self, dirc=None, file_eval=lambda a: True, dir_eval=lambda a: True):
        return local.ls(dirc, file_eval, dir_eval)
    
    def cp(self, f1, f2): copyfile(f1, f2)
   
    def __del__(self): pass

class Remote(remote.Shell, Base):
    
    def __init__(self, work_dir, ssh_config, *loaders):
        
        if ':' in work_dir:
            r = work_dir.split(':')
            server = r[0]
            dir_ = r[1]
            if dir_ == '': dir_ = '.'
        else:
            server = work_dir
            dir_ = '.'
            
        _gbl.logger.write('debug', f("Setting remote machine {r[0]} as host.."))
        
        remote.Shell.__init__(self, server, ssh_config)
        
        Base.__init__(self, dir_, *loaders)
    
    def hndl(self): return self.sftp
        
    def fileExists(self, file):
        try:
            self.sftp.stat(file)
            return True
        except FileNotFoundError:
           return False
    
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
    
    def open(self, name, s, prefetch=True):
        if 'r' in s: self.checkFile(name, throw=True)
        file = self.sftp.open(name, s)
        if 'w' not in s and prefetch: file.prefetch()
        return file
    
    def cp(self, f1, f2): self.run(f"cp {f1} {f2}")