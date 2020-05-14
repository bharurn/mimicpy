#rom ..utils.shell import Remote
#from collections import defaultdict
import re
import yaml
#from yaml.representer import Representer
import mimicpy._global as _global

class BaseCalc:
    
    def __init__(self, status=None):
        if not status:
            status = {'prepMM': '', 'prepQM': '', 'run': []}
        self._status = status
        self.log = ''
    
    def getcurrent(self, ext, level=False, exp=True):
        
        _dir = _global.host.pwd()+'/'
        
        if ext == 'top' or ext == 'mpt' or ext == 'mimic-tpr':
            return _dir+BaseCalc._getFile(self._status['prepMM'], ext)
        elif ext == 'mimic-tpr':
            return _dir+BaseCalc._getFile(self._status['prepQM'], ext)
        
        run = self._status['run']
        
        for i,d in enumerate(run[::-1]):
            ret = BaseCalc._getFile(d, ext)
            
            if ret != None:
                if level: return(i, _dir+ret)
                else: return _dir+ret
        
        if exp: raise Exception(f"Cannot find file type: {ext}")
    
    def __del__(self): self.saveToYaml()
        
    @staticmethod
    def _getFile(dirc, ext):
        lst = _global.host.ls(dirc=dirc, file_eval=lambda a: True if a.endswith(ext) else False, dir_eval=lambda a: False)
        
        if ext == 'cpt':
            lst = [l for l in lst if '_prev' not in l]
        
        if lst == []: return None
        elif len(lst) > 1: raise Exception("More than one file!")
        
        return dirc+'/'+lst[0]
    
    def getStatus(self): return self._status
    
    @classmethod
    def continueFrom(cls, session):
        if not hasattr(session, '_status'):
            raise Exception("Does not contain simulation status variables!")
        
        return cls(status=session._status)
    
    def saveToYaml(self):
        print(f"Saving status to _status.yaml..")
        y = yaml.dump(self._status)
        _global.host.write(y, '_status.yaml')
    
    @classmethod
    def continueFromYaml(cls):
       print(f"Loading session from _status.yaml..")
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
    
    #def moveMDResults(self, old, new):
    #    ls = self.ls(file_eval=lambda a: True if a.startswith(f"{old}.") or\
    #                 a.startswith(f"{old}_prev") else False, dir_eval = lambda a: False)
        
    #    for l in ls:
    #        a = l.split('.')
    #        n = a[0].replace(old, new)
            
    #        print(f"Renaming {l} to {n}.{a[-1]}")
    #        _global.host.cmd.rename(f"{l}", f"{n}.{a[-1]}")
    
    def dump(self, log_file):
        print(f"Dumping standard output from all MiMiC/MD runs so far to {log_file}..")
        
        _global.host.write(self.log, log_file)
    
    def _gmxhook(self, text):
        
        self.log += text
        
        BaseCalc._gmxerrhdnl(text)
        
        notes = BaseCalc._notes(text)
        
        if not notes.isspace():
            print(notes[:-1])
    
    @staticmethod
    def _gmxerrhdnl(text, dont_raise=False):
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
                err += f"Error running Gromacs!\n{err}"
            else:
                raise Exception(f"Error running Gromacs!\n{err}")
        
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
        
        return notes
    
    def gmx(self, cmd, **kwargs):
        
        if _global.gmx is None or _global.gmx.strip() == '':
            raise Exception('Gromacs executable not set!')
            
        gmx_ex = _global.gmx.strip()
        
        if 'nonverbose' not in kwargs:
            nonverbose = False
        else:
            nonverbose = True
            del kwargs['nonverbose']
            
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
        elif not nonverbose:
            print(f"Running {cmd}..")
        
        if mdrun: _global.host.runbg(cmd, hook=self._gmxhook, dirc=dirc)
        else: _global.host.run(cmd, stdin=stdin, hook=self._gmxhook, dirc=dirc)
    
    def cpmd(self, inp, out, onlycmd=False, noverbose=False, dirc=''):
        if _global.cpmd is None or _global.cpmd.strip() == '':
            raise Exception('CPMD executable not set!')
        
        if _global.cpmd_pp is None or _global.cpmd_pp.strip() == '':
            raise Exception('CPMD pseudopotential path not set!')
        
        cmd = f"{_global.cpmd} {inp} {_global.cpmd_pp} > {out}"
        
        if onlycmd: return cmd
        
        if not noverbose: print(cmd)
        
        _global.host.runbg(cmd, dirc=dirc)