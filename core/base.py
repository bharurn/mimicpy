from ..utils.shell import Remote
#from collections import defaultdict
import re
import yaml
#from yaml.representer import Representer
import mimicpy._global as _global

class BaseCalc:
    log = ''
    
    def __init__(self, status=[]):
        self._status = status
    
    def getcurrent(self, ext, level=False, exp=True):
        
        if ext == 'top': return self._status[0]+'/topol.top'
        elif ext == 'mpt': return self._status[0]+'/topol.mpt'
        
        for i,d in enumerate(self._status[::-1]):
            ret = BaseCalc._getFile(d, ext)
            
            if ret != None:
                if level: return(i, ret)
                else: return ret
        
        if exp: raise Exception(f"Cannot find file with extension {ext}")
    
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
    
    def saveToYaml(self, filename='_status.yaml'):
        #print(f"Saving session to {filename}..")
        #yaml.add_representer(defaultdict, Representer.represent_dict)
        _global.host.write(yaml.dump(self._status), filename)
    
    @classmethod
    def continueFromYaml(cls, filename='_status.yaml'):
       #print(f"Loading session from {filename}..")
       txt = _global.host.read(filename)
       _status = yaml.safe_load(txt)
       return cls(status=_status)
   
    def setcurrent(self, dirc=None):
        if dirc == None:
            dirc = self.dir
        _global.host.mkdir(dirc)
        if self._status[-1] != dirc:
            self._status.append(dirc)
    
    #def moveMDResults(self, old, new):
    #    ls = self.ls(file_eval=lambda a: True if a.startswith(f"{old}.") or\
    #                 a.startswith(f"{old}_prev") else False, dir_eval = lambda a: False)
        
    #    for l in ls:
    #        a = l.split('.')
    #        n = a[0].replace(old, new)
            
    #        print(f"Renaming {l} to {n}.{a[-1]}")
    #        _global.host.cmd.rename(f"{l}", f"{n}.{a[-1]}")
    
    @staticmethod
    def logToFile(self, log_file):
        print(f"Dumping standard output from all MiMiC/MD runs so far to {log_file}..")
        
        _global.host.write(BaseCalc.log, log_file)
    
    @staticmethod
    def _errorHandle(text, dont_raise=False):
        
        BaseCalc.log += text
        
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
    
    @staticmethod
    def gmx(cmd, **kwargs):
        
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
        
        cmd = f"{gmx_ex} {cmd}"
        
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
        
        text = _global.host.run(cmd, stdin=stdin, errorHandle=BaseCalc._errorHandle, query_rate = 3)
            
        notes = BaseCalc._notes(text)
        
        if notes != '' and not nonverbose:
            print('\n'+notes)
        
        return text
    
    @staticmethod
    def cpmd(inp, out, onlycmd=False, noverbose=False):
        if _global.cpmd is None or _global.cpmd.strip() == '':
            raise Exception('CPMD executable not set!')
        
        if _global.cpmd_pp is None or _global.cpmd_pp.strip() == '':
            raise Exception('CPMD pseudopotential path not set!')
        
        cmd = f"{_global.cpmd} {inp} {_global.cpmd_pp} > {out}"
        
        if onlycmd: return cmd
        
        if not noverbose: print(cmd)
        
        text = _global.host.run(cmd)
        
        return text