from .base import BaseCalc
import mimicpy._global as _global
from collections import defaultdict

class MD(BaseCalc):
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
                out = BaseCalc.gmx('mdrun', **kwargs, deffnm = new)
                return out
            else:
                return self.jobscript.add(BaseCalc.gmx('mdrun', **kwargs, deffnm = new, onlycmd=True))
        
        if 'noappend' in kwargs:
            if kwargs['noappend'] == True: out = self._do(new, **kwargs, noappend = '')
            else: out = _do(new, **kwargs)
        else:
            out = _do(new, **kwargs)
        
        if out == None:
            jid = _global.host.sbatch(self.jobscript)
            print("MDrun job submmitted as a Slurm job "
                 f"{self.jobscript.name}.sh with the job ID {jid}.."
                 f"\nThe host and/or this script can be safely closed..")
            return jid
        else:
            print("MDrun job submmitted on the current node..\n"
                 f"Please do not close host and/or this script until the job is done!!")
            return out
    
    def restart(self, new, until=0, extend=0, fromcrash=False, noappend=True):
        
        new_dir = '{new}/{new}'
        self.setcurrent(new)
        
        if fromcrash:
            return self.mdrun(new, s = self.getcurrent('tpr'), cpi = self.getcurrent('cpt'), noappend = noappend)
        elif until != 0:
            
            BaseCalc.gmx('convert-tpr', s = self.getcurrent('tpr'), until = until, o = f'{new_dir}.tpr')
            self.setcurrent('tpr', f"{new_dir}.tpr")
            return self.mdrun(s = f'{new_dir}.tpr', cpi = self.getcurrent('cpt'), deffnm = new_dir, noappend = noappend)
                
        elif extend:
            BaseCalc.gmx('convert-tpr', s = self.getcurrent('tpr'), extend = extend, o = f'{new_dir}.tpr')
            return self.mdrun(s = f'{new_dir}.tpr', cpi = '{cpt}.cpt', deffnm = new_dir, noappend = noappend)
        
        self.saveToYaml(new)
        
        
    def grompp(self, mdp, new, **kwargs):
        _global.host.write(str(mdp), f"{new}.mdp")
        
        gro = self.getcurrent('gro', level=True, exp=False)
        trr = self.getcurrent('trr', level=True, exp=False)
        
        if gro == None and trr ==  None:
            raise Exception("No coordinate data present..")
        
        if trr != None:
            if gro[0] < trr[0]:
                BaseCalc.gmx('grompp', f = f'{new}.mdp', c = gro[1],\
                     p = self.getcurrent('top'), o = f"{new}.tpr", **kwargs)
            else:
                BaseCalc.gmx('grompp', f = f'{new}.mdp', c = gro[1],\
                     p = self.getcurrent('top'), t=trr[1], o = f"{new}.tpr", **kwargs)
        else:
            BaseCalc.gmx('grompp', f = f'{new}.mdp', c = gro[1],\
                     p = self.getcurrent('top'), o = f"{new}.tpr", **kwargs)
        
    def run(self, mdp, **kwargs):
        new = mdp.name.tolower().replace(' ', '_')
        self.setdir(new)
        
        print(f"Starting classical MD calculation: {new}..")
        
        print(f"All files will be saved with the name {new}..")
        
        self.grompp(mdp, f'{new}/{new}', **kwargs)
        
        out = self.mdrun(f'{new}/{new}')
        
        self.saveToYaml(new)
        
        print("Done..")
        
        return out
        
    #def em(self, **kwargs):
    #    em_mdp = mdp.MDP.defaultEM()
    #   em_mdp.edit(**kwargs)
    #    return self.run(em_mdp)
    
    #def nvt(self, **kwargs):
    #    nvt_mdp = mdp.MDP.defaultNVT()
    #    nvt_mdp.edit(**kwargs)
    #    return self.run(nvt_mdp, r = self.getcurrent('coords'))
    
    #def npt(self, **kwargs):
    #    npt_mdp = mdp.MDP.defaultNPT()
    #    npt_mdp.edit(**kwargs)
    #    return self.run(npt_mdp, r = self.getcurrent('coords'))
    
    #def prd(self, **kwargs):
    #    prd_mdp = mdp.MDP.defaultPRD()
    #    prd_mdp.edit(**kwargs)
    #    return self.run(prd_mdp, r = self.getcurrent('coords'))
            
class MiMiC(MD):
    
    def setGMXSettings(self, **kwargs): self.gmx_opt = kwargs
    def setCPMDSettings(self, **kwargs): self.cpmd_opt = kwargs
    
    def run(self, inp, tpr=None):
        new = inp.info.lower().replace(' ', '_')
        self.setcurrent(new)
        _global.host.mkdir(f"{new}/cpmd")
        _global.host.mkdir(f"{new}/gmx")
        
        inp.mimic.paths = f"1\n{_global.host.pwd()+new}/gmx"
        _global.host.write(str(inp), f"{new}/cpmd/{new}.inp")
        
        if tpr is None: tpr = self.getcurrent('tpr')
        
        _global.host.cp(tpr, f'{new}/gmx/mimic.tpr')
        
        if self.jobscript:
            if not hasattr(self, 'gmx_opt'): self.gmx_opt = {}
            if not hasattr(self, 'cpmd_opt'): self.cpmd_opt = {}
            
            self.jobscript.add(self.gmx('mdrun', deffnm=f'{new}/gmx/mimic', onlycmd=True), **self.gmx_opt)
            self.jobscript.add(BaseCalc.cpmd(f"{new}/cpmd/{new}.inp", f"{new}/cpmd/{new}.out", onlycmd=True), **self.cpmd_opt)
            jid = _global.host.sbatch(self.jobscript)
            print("MiMiC run submmitted as a Slurm job "
                 f"{self.jobscript.name}.sh with the job ID {jid}.."
                 f"\nThe host and/or this script can be safely closed..")
            return jid
        else:
            out = self.gmx('mdrun', deffnm=f'{new}/gmx/mimic')
            BaseCalc.cpmd(f"{new}/cpmd/{new}.inp", f"{new}/cpmd/{new}.out")
            print("MiMiC run submmitted on the current node..\n"
                 f"Please do not close host and/or this script until the job is done!!")
            return out
        
        