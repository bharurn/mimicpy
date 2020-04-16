from pyshell.remote import ssh
import re

class GmxSSH(ssh.SSH):
    
    def __init__(self, server, gmx="gmx", directory="."):
        super().__init__(server)
        
        self.nongmxerr = True
        
        self.gmx = gmx
        self.log = ''
        ret = self.cd(directory, mkdir=True)
        
        if ret == -1:
            print("Setting current directory to $HOME..")
        else:
            if ret == 1:
                print(f"{directory} not found, creating new directory..")
            
            print(f"Setting current directory to {directory}..")
        
        self.query_rate = 3
    
    def _errorHandle(self, text, dont_raise=False):
        
        err = ''
        
        if 'ERROR' in text and self.nongmxerr:
            for line in text.splitlines():
                if 'ERROR' in line.upper():
                    if dont_raise:
                        err += line
                    else:
                        raise Exception(line)
        
        if 'Halting program' in text or 'Fatal error' in text:
            
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
        
    def _notes(self, text):
        
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
     
    def run(self, cmd, stdin=None, parallel=False, errorHandle=None, onNewChan=False):
        splt = cmd.split()
        
        if 'gmx' in splt[0]:    
            print(f"Running Gromacs {splt[1]}..")
            self.query_rate = 3
            self.nongmxerr = False
        else:
            self.nongmxerr = True
            self.query_rate = 0.3
        
        if errorHandle == None:
            errorHandle = self._errorHandle
            
        text = super().run(cmd, stdin=stdin, parallel=parallel, errorHandle=errorHandle, onNewChan=onNewChan)
        
        notes = self._notes(text)
            
        if 'gmx' in splt[0]:    
            self.log += f"========================\n{cmd}\n========================\n"+text+'\n\n'
            if notes != '':
                print('\n'+notes)
        
        return text
            
    def restartMD(self, old, new, job, until=0, extend=0, fromcrash=False):
        
        if isinstance(old, str):
            tpr = old
            cpt = old
        else:
            tpr = old[0]
            cpt = old[1]
        
        if fromcrash:
            job.add(f'{self.gmx} mdrun -s {tpr}.tpr -cpi {cpt}.cpt -noappend -deffnm {new}')
        elif until:
            job.add(f'{self.gmx} convert-tpr -s {tpr}.tpr -until {until} -o {new}.tpr')
            job.add(f'{self.gmx} mdrun -s {new}.tpr -cpi {cpt}.cpt -noappend -deffnm {new}')
        elif extend:
            job.add(f'{self.gmx} convert-tpr -s {tpr}.tpr -extedn {extend} -o {new}.tpr')
            job.add(f'{self.gmx} mdrun -s {new}.tpr -cpi {cpt}.cpt -noappend -deffnm {new}')
        
        return self.sbatch(job)
    
    def moveMDResults(self, old, new):
        ls = self.ls(file_eval=lambda a: True if a.startswith(old) else False, dir_eval = lambda a: False)
        
        for l in ls:
            a = l.split('.')
            n = a[0].replace(old, new)
            
            print(f"Renaming {l} to {n}.{a[-1]}")
            self.sftp.rename(f"{l}", f"{n}.{a[-1]}")