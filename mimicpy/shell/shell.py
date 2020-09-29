from . import _local as local
import os
from stat import S_ISDIR, S_ISREG
from shutil import copyfile
import re
from .._global import _Global as _gbl
from abc import ABC, abstractmethod

class Base(ABC):

    @abstractmethod
    def __init__(self, directory, *loaders):

        if directory.strip() == '':
            directory = '.'

        ret = self.cd(directory, mkdir=True)

        if ret == -1: # dir is home directory
            pass
        else:
            if ret == 1: # dir doesn't exist and was created by cd()
                _gbl.logger.write('debug', f"{directory} not found, creating new directory..")
            # if ret == 0: --> dir already exists
            _gbl.logger.write('debug', f"Setting current working directory to {self.pwd()}..")

        self.loaders = [] # list of loaders
        self.loader_str = '' # stringified version of loaders list
        self.loader_out = '' # output of running loaders

        if loaders: # add loaders to self.loaders
            self.addLoaders(*loaders)

    def addLoaders(self, *loader):
         self.loaders.extend(loader) # add loaders
         if self.loader_str.strip() == '': # check if loaders already there
             self.loader_str += ' ; '.join(loader) # join all loaders with ';'
         else:
             self.loader_str = ' ; '.join(loader) # join all loaders with ';'

         # replace defaultHook in run(), so it doesn't check for errors
         print_hook = lambda cmd, out: _gbl.logger.write('debug2', f'Running loader {cmd}..\n{out}')

         # get output of loaders, so we can remove it from all subsequent run() outputs
         self.loader_out = self.run(self.loader_str, hook=print_hook, fresh=True)

    #convenience functions
    def rename(self, a, b): self.hndl().rename(a, b)
    def rm(self, a): self.hndl().remove(a)
    def pwd(self): return self.hndl().getcwd()
    def isLocal(self):
        if self.name == 'localhost': return True
        else: return False

    @abstractmethod
    def hndl(self): pass

    def checkFile(self, file, throw=False):
        """ Check if file exists, for debugging """
        ret = self.fileExists(file)

        if ret == True:
            pass
            #_gbl.logger.write('debug', f"{file} found!")
        elif throw:
            raise FileNotFoundError(f"{file} not found!")

    def read(self, file, asbytes=False):
        """ Read and return contents of full file """
        mode = 'rb' if asbytes else 'r'
        with self.open(file, mode) as f:
            out = f.read()
            if not asbytes:
                # sftp file objects, when read, returns only bytes
                # so need to decode it if self is remote
                try:
                    return out.decode()
                except (UnicodeDecodeError, AttributeError): # if error, means out is string
                    return out
            else:
                # mode should be rb
                # so o/p will be in bytes, just return
                return out

    def write(self, content, file, asbytes=False):
        """ Write file """
        mode = 'wb' if asbytes else 'w'
        with self.open(file, mode) as f: f.write(content)

    def mkdir(self, directory):
        """ Make directory """
        if not self.fileExists(directory):
            self.hndl().mkdir(directory)
            return 0
        else: return 1

    def cd(self, directory, mkdir=False):
        """
        Change directory
        If mkdir is True, cd() will make a new directory if not found
        Otherwise it will raise exception
        """
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

    def sbatch(self, job, dirc=''):
        """
        Write jobscript to file
        And run sbatch <jobscript.sh>
        """
        jbs = f'{dirc}/{job.name}.sh'

        self.write(str(job), jbs)

        if job.noCommands():
            pass
#            raise SlurmBatchError(jbs, "No commands found!")

        jid = 0
        def _sbatch_err(txt):
            if 'error' in txt.lower():
                pass
#                raise SlurmBatchError(jbs, txt)
            else:
                match = re.search(r'Submitted batch job (\w*)', txt)
                if match:
                    global jid
                    jid = match.groups()[0]
                else:
                    pass
#                    raise SlurmBatchError(jbs, txt)

        self.run(f'sbatch {job.name}.sh', hook=_sbatch_err, dirc=dirc, fresh=True)

        return jid

    def scancel(self, jobid):
        """Cancel slurm job"""
        return self.run(f'scancel {jobid}', fresh=True)

class Local(Base):
    """

    Class to handle localhost specific shell function

    """

    def __init__(self, directory, shell_path, *loaders):
        """Init cwd and shell executaable path"""
        self.name = 'localhost'
        if shell_path:
            self.shell_path = shell_path
        else: # by default uses $SHELL to get shell path
            self.shell_path = os.environ['SHELL'] # check env vars and get shell path, works only for UNIX
        super().__init__(directory, *loaders)

    # return handle for file io
    def hndl(self): return os

    def fileExists(self, file):
        """Check if file exists"""
        return os.path.isfile(file) or os.path.isdir(file)

    def open(self, file, mode):
        """Wrapper around open()"""
        if 'r' in mode: self.checkFile(file, throw=True)
        return open(file, mode)

    def run(self, cmd, stdin=None, hook=None, fresh=False, dirc=''):
        """Local shell specific run command"""
        replace = ''
        if not fresh and self.loader_str:
            cmd = self.loader_str + ' ; ' + cmd # add loader string
            replace = self.loader_out # remove loader_out cmd from output
        if dirc != '': cmd = f'cd {dirc} ; ' + cmd # move to dirc if specified
        return local.run(cmd, self.shell_path, replace, stdin=stdin, hook=hook) # call run from _local.py

    def runbg(self, cmd, hook=None, fresh=False, dirc='', query_rate=1):
        """Local shell specific run background command"""
        replace = ''
        if not fresh:
            cmd = self.loader_str + ' ; ' + cmd
            replace = self.loader_out
        if not dirc != '': cmd = f'cd {dirc} ; ' + cmd
        return local.runbg(cmd, self.shell_path, replace, hook=hook) # call runbg from _local.py

    def ls(self, dirc=None, file_eval=lambda a: True, dir_eval=lambda a: True):
        """Local shell specific ls command to return files/folders as list"""
        return local.ls(dirc, file_eval, dir_eval)

    # copy command
    def cp(self, f1, f2): copyfile(f1, f2)

    def join(self, *args): return os.path.join(*args)

    def dirname(self, file): return os.path.dirname(file)

    # provided to match remote.__del__()
    def __del__(self): pass
