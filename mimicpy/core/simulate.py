from .base import BaseHandle
from .._global import _Global as _global
from collections import defaultdict
from ..utils.errors import MiMiCPyError

class MD(BaseHandle):
    def __init__(self, status=defaultdict(list), settings=None):
        super().__init__(status)
        self.jobscript = None
        if settings: self.setSlurmSettings()
    
    def setSlurmSettings(self, settings):
        _global.logger.write('debug', f"Setting Slurm job settings from jobscript {settings.name}..")
        self.jobscript = settings
        if _global.host.loaders != []:
            _global.logger.write('debug', f"Transferring loader commands from host to job script..")
            self.jobscript.addMany(_global.host.loaders)
        
        return self.jobscript
            
    def mdrun(self, new, **kwargs):
        
        def _do(new, **kwargs):
            if self.jobscript is None:
                self.gmx('mdrun', **kwargs, deffnm = new)
            else:
                self.jobscript.add(self.gmx('mdrun', **kwargs, deffnm = new, onlycmd=True))
        
        if 'noappend' in kwargs:
            if kwargs['noappend'] == True:
                del kwargs['noappend']
                _do(new, **kwargs, noappend = '')
            else: _do(new, **kwargs)
        else:
            _do(new, **kwargs)
        
        if self.jobscript:
            if 'dirc' in kwargs:
                dirc = kwargs['dirc']
                del kwargs['dirc']
            else:
                dirc = ''
                    
            jid = _global.host.sbatch(self.jobscript, dirc=dirc)
            _global.logger.write('info', "Gromacs simulations submmitted as a Slurm job "
                 f"{self.jobscript.name}.sh with the job ID {jid}.."
                 f"\nThe host and/or this script can be safely closed..")
            return jid
        else:
            _global.logger.write('info', "Gromacs simulation is now running as a background job..")
            if _global.host.name != 'localhost':
                _global.logger.write('info', "Please do not close remote host until the job is done!")
    
    def restart(self, new, until=0, extend=0, fromcrash=False, noappend=True):
        
        self.setcurrent(new)
        
        if fromcrash:
            out = self.mdrun(new, s = self.getcurrent('tpr'), cpi = self.getcurrent('cpt'), noappend = noappend, dirc=new)
        elif until != 0:
            
            self.gmx('convert-tpr', s = self.getcurrent('tpr'), until = until, o = f'{new}.tpr', dirc=new)
            out = self.mdrun(new, s = f'{new}.tpr', cpi = self.getcurrent('cpt'), noappend = noappend, dirc=new)
                
        elif extend:
            self.gmx('convert-tpr', s = self.getcurrent('tpr'), extend = extend, o = f'{new}.tpr', dirc=new)
            out = self.mdrun(new, s = f'{new}.tpr', cpi = '{cpt}.cpt', noappend = noappend, dirc=new)
        
        self.saveToYaml()
        return out
        
        
    def grompp(self, mdp, new, **kwargs):
        new = new.lower().replace(' ', '_')
        if 'dirc' in kwargs:
            mdp_file = f'{kwargs["dirc"]}/{new}.mdp'
        else:
            mdp_file = f'{new}.mdp'
        _global.host.write(str(mdp), mdp_file)
        
        if 'gro' in kwargs:
            gro_file = kwargs['gro']
            del kwargs['gro']
        elif 'trr' in kwargs:
            trr_file = kwargs['trr']
            del kwargs['trr']
        
        gro = self.getcurrentNone(gro_file, 'gro', level=True, exp=False)
        trr = self.getcurrentNone(trr_file, 'trr', level=True, exp=False)
        
        if gro == None and trr ==  None:
            raise MiMiCPyError("No coordinate data (gro/trr file) was found..")
        
        if trr != None:
            if gro[0] < trr[0]:
                self.gmx('grompp', f = f'{new}.mdp', c = gro[1],\
                     p = self.getcurrent('top'), o = f"{new}.tpr", **kwargs)
            else:
                self.gmx('grompp', f = f'{new}.mdp', c = gro[1],\
                     p = self.getcurrent('top'), t=trr[1], o = f"{new}.tpr", **kwargs)
        else:
            self.gmx('grompp', f = f'{new}.mdp', c = gro[1],\
                     p = self.getcurrent('top'), o = f"{new}.tpr", **kwargs)
        
    def run(self, mdp, **kwargs):
        new = mdp.name.lower().replace(' ', '_')
        
        if 'dir' in kwargs: # dir is an argument for the directory, if dir = None, then use cwd
            if kwargs['dir'] == None: kwargs['dir'] = ''
            _dir = kwargs['dir']
        else:   
            _dir = new
        self.setcurrent(_dir)
        
        _global.logger.write('info', f"Starting classical MD calculation: {new}..")
        
        _global.logger.write('debug', f"All files will be saved to the directory {_dir}..")
        
        self.grompp(mdp, f'{new}', dirc=_dir, **kwargs)
        
        out = self.mdrun(f'{new}', dirc=_dir)
        
        self.saveToYaml()
        
        return out
        
class MiMiC(MD):
    
    def setGMXSettings(self, **kwargs): self.gmx_opt = kwargs
    def setCPMDSettings(self, **kwargs): self.cpmd_opt = kwargs
    
    def run(self, inp, tpr=None, dirc=''):
        new = inp.info.lower().replace(' ', '_')
        self.setcurrent(new)
        _global.host.mkdir(f"{new}/cpmd")
        _global.host.mkdir(f"{new}/gmx")
        
        inp.mimic.paths = f"1\n{_global.host.pwd()+new}/gmx"
        _global.logger.write('debug2', "Set PATH in MIMIC section as {inp.mimic.paths}")
        _global.host.write(str(inp), f"{new}/cpmd/{new}.inp")
        
        tpr = self.getcurrentNone('mimic-tpr')
        
        _global.host.cp(tpr, f'{new}/gmx/mimic.tpr')
        
        if self.jobscript:
            if not hasattr(self, 'gmx_opt'): self.gmx_opt = {}
            if not hasattr(self, 'cpmd_opt'): self.cpmd_opt = {}
            
            self.jobscript.add(self.gmx('mdrun', deffnm=f'gmx/mimic', onlycmd=True, dirc='new'), **self.gmx_opt)
            self.jobscript.add(self.cpmd(f"cpmd/{new}.inp", f"cpmd/{new}.out", onlycmd=True, dirc='new'), **self.cpmd_opt)
            jid = _global.host.sbatch(self.jobscript, dirc=dirc)
            _global.logger.write('info', "MiMiC run submmitted as a Slurm job "
                 f"{self.jobscript.name}.sh with the job ID {jid}.."
                 f"\nThe host and/or this script can be safely closed..")
            return jid
        else:
            self.gmx('mdrun', deffnm=f'gmx/mimic', dirc='new')
            self.cpmd(f"cpmd/{new}.inp", f"cpmd/{new}.out", dirc='new')
            _global.logger.write('info', "MiMiC simulation is now running in the background..")
            
            if _global.host.name != 'localhost':
                _global.logger.write('info', "Please do not close remote host until the job is done!")
            
        self.saveToYaml()
        
        