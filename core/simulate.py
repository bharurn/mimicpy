from ..utils.scripts import mdp
from .base import Run
import mimicpy._global as _global
from collections import defaultdict

class MD(Run):
    def __init__(self, status=defaultdict(list), settings=None):
        super().__init__(status)
        if settings: self.setSlurmSettings()
    
    def setSlurmSettings(self, settings):
        print(f"Setting Slurm job settings from jobscript {settings.name}..")
        self.jobscript = settings
        print(f"Transferring sources, modules and other header commands from host to job script..")
        self.jobscript.module(_global.host.modules)
        self.jobscript.source(_global.host.sources)
        self.jobscript.addMany(_global.host.loaders)
        
        return self.jobscript
            
    def mdrun(self, new, **kwargs):
        
        def _do(new, **kwargs):
            if self.jobscript is None:
                out = Run.gmx('mdrun', **kwargs, deffnm = new)
                return out
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
        self.saveToYaml()
        
        if out == None:
            jid = _global.host.sbatch(self.jobscript)
            print("MDrun job submmitted as a Slurm job "
                 f"{self.jobscript.name}.sh with the job ID {jid}.."
                 f"\nThe host and/or this script can be safely closed..")
            return jid
        else:
            print("MDrun job submmitted on the current login node"
                 f"Please do not close this script or config.host until the job is done!!")
            return out
    
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
        
        print(f"Using {self.getcurrent('coords')}, {self.getcurrent('topology')} for grompp")
        
        Run.gmx('grompp', f = f'{new}.mdp', c = self.getcurrent('coords'),\
                     p = self.getcurrent('topology'), o = f"{new}.tpr", **kwargs)
        
        self.setcurrent('tpr', f"{new}.tpr")
        
        out = self.mdrun(new)
        
        print("Done..")
        
        return out
        
    def em(self, **kwargs):
        em_mdp = mdp.MDP.defaultEM()
        em_mdp.edit(**kwargs)
        return self.calc(em_mdp, 'em', disp='Minimization')
    
    def nvt(self, **kwargs):
        nvt_mdp = mdp.MDP.defaultNVT()
        nvt_mdp.edit(**kwargs)
        return self.calc(nvt_mdp, 'nvt', disp='NVT simulation', r = self.getcurrent('coords'))
    
    def npt(self, **kwargs):
        npt_mdp = mdp.MDP.defaultNPT()
        npt_mdp.edit(**kwargs)
        return self.calc(npt_mdp, 'npt', disp='NPT simulation', r = self.getcurrent('coords'))
    
    def prd(self, **kwargs):
        prd_mdp = mdp.MDP.defaultPRD()
        prd_mdp.edit(**kwargs)
        return self.calc(prd_mdp, 'prd', disp='Production run', r = self.getcurrent('coords'))
            