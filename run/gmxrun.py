from .. import host
from ..shell import remote
from collections import defaultdict
import re

class GMX:
    log = ''
    def __init__(self):
        self._status = defaultdict(list)
    
    def setcurrent(self, key, val): self._status[key].append(val)
    def getcurrent(self, key): return self._status[key][-1]
    def gethistory(self, key): return self._status[key]
    
    def continueFrom(self, session):
        if hasattr(session, '_status'):
            raise Exception("Does not contain simulation status variables!")
        
        self._status = session._status
    
    def saveSessionToFile(self, filename):
        pass
    
    @classmethod
    def continueSessionFromFile(cls, filename):
        pass
    
    def moveMDResults(self, old, new):
        ls = self.ls(file_eval=lambda a: True if a.startswith(f"{old}.") or\
                     a.startswith(f"{old}_prev") else False, dir_eval = lambda a: False)
        
        for l in ls:
            a = l.split('.')
            n = a[0].replace(old, new)
            
            print(f"Renaming {l} to {n}.{a[-1]}")
            host.cmd.rename(f"{l}", f"{n}.{a[-1]}")
    
    def logToFile(self, log):
        print(f"Dumping standard output from all Gromacs runs so far to {log}..")
        
        host.cmd.write(host.cmd.log, log)
    
    @staticmethod
    def _errorHandle(text, dont_raise=False):
        
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
        
        if not hasattr(host, 'gmx'):
            raise Exception('Gromacs executable not set! Please set it in host..')
            
        gmx_ex = host.gmx.strip()
        
        if gmx_ex[:3] != 'gmx':
            raise Exception(f'{host.gmx} is an invalid Gromacs executable! Please set it correctly in host..')
        
        splt = cmd.split()
        
        if 'noverbose' not in kwargs:
            print(f"Running Gromacs {splt[0]}..")
        else:
            del kwargs['noverbose']
            
        if type(host.cmd) is remote.SSH: host.cmd.query_rate = 3
        
        cmd = f"{gmx_ex} {cmd}"
        
        stdin = None
        onlycmd = False
        
        for k, v in kwargs.items():
            if k != 'stdin':
                cmd += f" -{k} {v}"
            elif k == 'stdin':
                stdin = v
            elif k == 'onlycmd':
                onlycmd = True
                
        if onlycmd: return cmd
        
        text = host.cmd.run(cmd, stdin=stdin, errorHandle=GMX._errorHandle)
        
        notes = GMX._notes(text)
            
        GMX.log += f"========================\n{cmd}\n========================\n"+text+'\n\n'
        if notes != '':
            print('\n'+notes)
        
        return text