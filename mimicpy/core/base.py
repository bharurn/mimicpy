import re
import yaml
from .._global import _Global as _global
from ..utils.errors import GromacsError, MiMiCPyError, EnvNotSetError, defaultHook
from ..utils.logger import Logger, LogString
import sys

class BaseHandle:
    
    def __init__(self, status=None):
        if not status:
            status = {'prepMM': '', 'prepQM': '', 'run': ['']}
        self._status = status
        self.log = LogString()
        self.logger = Logger(log=_global.host.open('gmx.log', 'w'), notes=sys.stdout)
        self.current_cmd = 'gmx'
    
    def getcurrent(self, ext, level=False, exp=True):
        
        _dir = _global.host.pwd()+'/'
        
        if ext == 'top' or ext == 'mpt' or ext == 'mimic-tpr':
            return _dir+BaseHandle._getFile(self._status['prepMM'], ext)
        elif ext == 'mimic-tpr':
            return _dir+BaseHandle._getFile(self._status['prepQM'], ext)
        
        run =  [self._status['prepMM'], self._status['prepQM']] + self._status['run']
        
        for i,d in enumerate(run[::-1]):
            ret = BaseHandle._getFile(d, ext)
            
            if ret != None:
                if level: return(i, _dir+ret)
                else: return _dir+ret
        
        if exp: raise FileNotFoundError(f"Cannot find file with extension {ext}")
    
    def __del__(self):
        self.logger.close()
        self.toYaml()
        
    @staticmethod
    def _getFile(dirc, ext):
        lst = _global.host.ls(dirc=dirc, file_eval=lambda a: True if a.endswith(ext) else False, dir_eval=lambda a: False)
        
        if ext == 'cpt':
            lst = [l for l in lst if '_prev' or '_step' not in l]
        
        if lst == []: return None
        elif len(lst) > 1: raise MiMiCPyError(f"More than one current file found with extension {ext}!"
                                        f"\nFiles found: {','.join(lst)}")
        if dirc.strip() != '': dirc += '/'
            
        return dirc+lst[0]
    
    def getStatus(self): return self._status
    
    @classmethod
    def continueFrom(cls, session):
        if not hasattr(session, '_status'):
            raise MiMiCPyError("No simulation status variables were found for session!")
        
        return cls(status=session._status)
    
    def toYaml(self):
        _global.logger.write('debug', f"Saving status to _status.yaml..")
        y = yaml.dump(self._status)
        _global.host.write(y, '_status.yaml')
    
    @classmethod
    def fromYaml(cls):
       _global.logger.write('debug', f"Loading status from _status.yaml..")
       txt = _global.host.read('_status.yaml')
       _status = yaml.safe_load(txt)
       return cls(status=_status)
   
    def setcurrent(self, dirc=None, key='run'):
        if dirc == None:
            dirc = self.dir
        _global.host.mkdir(dirc)
        if key == 'prepMM' or key == 'prepQM':
            self._status[key] = dirc
        elif dirc not in self._status['run']:
            self._status['run'].append(dirc)
    
    def _gmxhook(self, cmd, text):
        defaultHook(cmd, text)
        
        self.logger.write('log', f"============Running {cmd}============\n")
        self.logger.write('log', text)
        
        self._gmxerrhdnl(cmd, text)
        
        notes = BaseHandle._notes(text)
        
        if not notes.isspace() and notes.strip() != '':
            self.logger.write('notes', f"From {self.current_cmd}")
            self.logger.write('notes', notes)
    
    @staticmethod
    def _gmxerrhdnl(gmx_cmd, text, dont_raise=False):
        if 'Halting program' in text or 'Fatal error' in text or 'Error in user input' in text:
            
            pattern = re.compile("^-+$")
            
            msg = []
            start = False
            
            for line in text.splitlines()[::-1]:
                if line.startswith('For more information and tips for troubleshooting') and not start:
                    start=True
                elif start:
                    if pattern.match(line):
                        break
                    else:
                        msg.append(line)
                    
                err = '\n'.join(msg[::-1])
                  
            if dont_raise:
                err += f"Error running Gromacs {gmx_cmd}!\n{err}"
            else:
                raise GromacsError(gmx_cmd, err)
        
        if dont_raise:
            return err
        else:
            return None
        
    @staticmethod
    def _notes(text):
        
        pattern = re.compile("^-+$")
        
        notes = ''
        start = False
        for i, line in enumerate(text.splitlines()):
            if line.startswith('NOTE') or line.startswith('WARNING') or 'One or more water molecules can not be settled.' in line:
                notes += line + '\n'
                start = True
            elif start:
                if line.strip() == '' or line.startswith('WARNING') or line.startswith('NOTE') or pattern.match(line):
                    start = False
                else:
                    notes += line + '\n'
        
        return ''.join(filter(lambda x: not re.match(r'^\s*$', x), notes)) # remove blank lines
    
     
    def gmx(self, cmd, **kwargs):
        
        if _global.gmx is None or _global.gmx.strip() == '':
            raise EnvNotSetError('Gromacs executable', 'gmx')
            
        gmx_ex = _global.gmx.strip()
        
        if 'onlycmd' in kwargs:
            onlycmd = True
            del kwargs['onlycmd']
        else:
            onlycmd = False
        
        if 'dirc' in kwargs:
            dirc = kwargs['dirc']
            del kwargs['dirc']
        else:
            dirc = ''
        
        if cmd == 'mdrun': mdrun = True
        else: mdrun = False
        
        cmd = f"{gmx_ex} -quiet {cmd}"
        
        stdin = None
        
        for k, v in kwargs.items():
            if k != 'stdin':
                cmd += f" -{k} {v}"
            elif k == 'stdin':
                stdin = v
                
        if onlycmd:
            return cmd
        
        _global.logger.write('debug', f"Running {cmd}..")
        
        self.current_gmx = cmd
        if mdrun: _global.host.runbg(cmd, hook=self._gmxhook, dirc=dirc)
        else: _global.host.run(cmd, stdin=stdin, hook=self._gmxhook, dirc=dirc)
    
    def cpmd(self, inp, out, onlycmd=False, noverbose=False, dirc=''):
        if _global.cpmd is None or _global.cpmd.strip() == '':
            raise EnvNotSetError('CPMD executable', 'cpmd')
        
        if _global.cpmd_pp is None or _global.cpmd_pp.strip() == '':
            raise EnvNotSetError('CPMD pseudopotential', 'cpmd_pp')
        
        cmd = f"{_global.cpmd} {inp} {_global.cpmd_pp} > {out}"
        
        if onlycmd: return cmd
        
        if not noverbose: print(cmd)
        
        _global.logger.write('debug', "Running {cmd}..")
        _global.host.runbg(cmd, dirc=dirc)