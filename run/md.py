from pygmx.run import mdp, gmxrun
from pygmx import host

class MD(gmxrun.GMX):
    
    def __init__(self, **kwargs):
        
        self.jobscript = None
        
        for k,v in kwargs.items():
            self.setcurrent(k, v)
    
    def SlurmSettings(self, settings):
        print(f"Setting SLURM job settings from jobscript {settings.name}..")
        self.jobscript = settings
        print(f"Transferring sources, modules and other header commands from Host to job script..")
        self.jobscript.module(host.cmd.modules)
        self.jobscript.source(host.cmd.sources)
        self.jobscript.addMany(host.cmd.loaders)
        
        return self.jobscript
            
    def mdrun(self, new, **kwargs):
        
        def _do(self, new, **kwargs):
            if self.jobscript in None:
                self.gmx('mdrun', **kwargs, deffnm = new)
            else:
                self.jobscript.add(self.gmx('mdrun', **kwargs, deffnm = new, onlycmd=True))
        
        if kwargs['noappend'] == True:
            out = self._do(new, **kwargs, noappend = '')
        else:
            out = self._do(new, **kwargs)
            
        self.setcurrent('edr', f"{new}.edr")
        self.setcurrent('trr', f"{new}.trr")
        self.setcurrent('cpt', f"{new}.cpt")
        self.setcurrent('log', f"{new}.log")
        
        if out == None:
            print("MDrun job submmitted as a SLURM job script"
                 f"{self.jobscript.name}.sh..\nPlease use squeue() for details and check for job ID {out}..")
            return host.cmd.sbatch(self.jobscript)
        else:
            print("Running mdrun on the login node..\nConsider using SLURM if this is a production run..")
            return None
    
    def continueRun(self, new, until=0, extend=0, fromcrash=False, noappend=True):
        
        if fromcrash:
            return self.mdrun(new, s = self.getcurrent('tpr'), cpi = self.getcurrent('cpt'), noappend = noappend)
        elif until != 0:
            
            self.gmx('convert-tpr', s = self.getcurrent('tpr'), until = until, o = f'{new}.tpr')
            self.setcurrent('tpr', f"{new}.tpr")
            return self.mdrun(s = f'{new}.tpr', cpi = self.getcurrent('cpt'), deffnm = new, noappend = noappend)
                
        elif extend:
            self.gmx('convert-tpr', s = self.getcurrent('tpr'), extend = extend, o = f'{new}.tpr')
            self.setcurrent('tpr', f"{new}.tpr")
            return self.mdrun(s = f'{new}.tpr', cpi = '{cpt}.cpt', deffnm = new, noappend = noappend)
    
    def calc(self, mdp, new, **kwargs):
        
        if 'disp' not in kwargs: disp = new.title()
        else:
            disp = kwargs['name']
            del kwargs['name']
        
        print(f"**Starting classical MD calculation: {disp}**")
        
        print(f"All files will be saved with the name {new}..")
        
        self.write(str(mdp), f"{new}.mdp")
        
        print("Running grompp..")
        
        self.gmx('grompp', f = f'{new}.mdp', c = self.getcurrent('coords'),\
                     p = self.getcurrent('topology'), o = f"{new}.tpr", **kwargs)
        
        self.setcurrent('tpr', f"{new}.tpr")
        
        out = self.mdrun(new)
        
        print("**Done**")
        
        return out
        
    def em(self, em_mdp = mdp.MDP.defaultEM()): return self.calc(em_mdp, 'em', disp='Minimization')
    
    def nvt(self, nvt_mdp = mdp.MDP.defaultNVT()): return self.calc(nvt_mdp, 'nvt', disp='NVT simulation', r = self.current('coords'))
    
    def npt(self, npt_mdp = mdp.MDP.defaultNPT()): return self.calc(npt_mdp, 'npt', disp='NPT simulation', r = self.current('coords'))
    
    def prd(self, prd_mdp = mdp.MDP.defaultPRD()): return self.calc(prd_mdp, 'prd', disp='Production run', r = self.current('coords'))
            