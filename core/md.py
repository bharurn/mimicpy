from ..utils.scripts import mdp
from .base import Run
import mimicpy._global as _global

class MD(Run):
    
    def SlurmSettings(self, settings):
        print(f"Setting Slurm job settings from jobscript {settings.name}..")
        self.jobscript = settings
        print(f"Transferring sources, modules and other header commands from _global.host to job script..")
        self.jobscript.module(_global.host.modules)
        self.jobscript.source(_global.host.sources)
        self.jobscript.addMany(_global.host.loaders)
        
        return self.jobscript
            
    def mdrun(self, new, **kwargs):
        
        def _do(new, **kwargs):
            if self.jobscript is None:
                 # shell.local -- add polling for stdout and print timestep
                _global.host.redirectStdout(f"{new}.log")
                Run.gmx('mdrun', **kwargs, deffnm = new)
                print("Done..")
            else:
                return self.jobscript.add(Run.gmx('mdrun', **kwargs, deffnm = new, onlycmd=True, nonverbose=True))
        
        if 'noappend' in kwargs:
            if kwargs['noappend'] == True: out = self._do(new, **kwargs, noappend = '')
            else: out = _do(new, **kwargs)
        else:
            out = _do(new, **kwargs)
            
        self.setcurrent('edr', f"{new}.edr")
        self.setcurrent('trr', f"{new}.trr")
        self.setcurrent('cpt', f"{new}.cpt")
        self.setcurrent('log', f"{new}.log")
        
        if out == None:
            print("MDrun job submmitted as a Slurm job script"
                 f"{self.jobscript.name}.sh with the job ID {out}.."
                 f"\nThe host and/or this script can be safely closed..")
            return _global.host.sbatch(self.jobscript)
        else:
            return None
    
    def continueRun(self, new, until=0, extend=0, fromcrash=False, noappend=True):
        
        if fromcrash:
            return self.mdrun(new, s = self.getcurrent('tpr'), cpi = self.getcurrent('cpt'), noappend = noappend)
        elif until != 0:
            
            Run.gmx('convert-tpr', s = self.getcurrent('tpr'), until = until, o = f'{new}.tpr')
            self.setcurrent('tpr', f"{new}.tpr")
            return self.mdrun(s = f'{new}.tpr', cpi = self.getcurrent('cpt'), deffnm = new, noappend = noappend)
                
        elif extend:
            Run.gmx('convert-tpr', s = self.getcurrent('tpr'), extend = extend, o = f'{new}.tpr')
            self.setcurrent('tpr', f"{new}.tpr")
            return self.mdrun(s = f'{new}.tpr', cpi = '{cpt}.cpt', deffnm = new, noappend = noappend)
    
    def calc(self, mdp, new, **kwargs):
        
        if 'disp' not in kwargs: disp = new.title()
        else:
            disp = kwargs['disp']
            del kwargs['disp']
        
        print(f"Starting classical MD calculation: {disp}..")
        
        print(f"All files will be saved with the name {new}..")
        
        _global.host.write(str(mdp), f"{new}.mdp")
        
        Run.gmx('grompp', f = f'{new}.mdp', c = self.getcurrent('coords'),\
                     p = self.getcurrent('topology'), o = f"{new}.tpr", **kwargs)
        
        self.setcurrent('tpr', f"{new}.tpr")
        
        out = self.mdrun(new)
        
        print("Done..")
        
        return out
        
    def em(self, em_mdp = mdp.MDP.defaultEM()): return self.calc(em_mdp, 'em', disp='Minimization')
    
    def nvt(self, nvt_mdp = mdp.MDP.defaultNVT()): return self.calc(nvt_mdp, 'nvt', disp='NVT simulation', r = self.current('coords'))
    
    def npt(self, npt_mdp = mdp.MDP.defaultNPT()): return self.calc(npt_mdp, 'npt', disp='NPT simulation', r = self.current('coords'))
    
    def prd(self, prd_mdp = mdp.MDP.defaultPRD()): return self.calc(prd_mdp, 'prd', disp='Production run', r = self.current('coords'))
            