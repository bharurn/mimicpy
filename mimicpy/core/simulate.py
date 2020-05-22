# This is part of MiMiCPy

"""

This module contains handles to run MD and MiMiC simulations

"""


from .base import BaseHandle
from .._global import _Global as _global
from collections import defaultdict
from ..utils.errors import MiMiCPyError

class MD(BaseHandle):
    """
    Simulte MD runs using grompp, convetr-tpr, mdrun
    Inherits from .core.base.BaseHandle
    
    """
    def __init__(self, status=defaultdict(list), settings=None):
        """Class constructor"""
        super().__init__(status)
        self.jobscript = None
        if settings: self.setSlurmSettings()
    
    def setSlurmSettings(self, settings):
        """Set the Slurm settings from a jobscript"""
        _global.logger.write('debug', f"Setting Slurm job settings from jobscript {settings.name}..")
        self.jobscript = settings
        if _global.host.loaders != []: # transfer host loaders to jobscript
            _global.logger.write('debug', f"Transferring loader commands from host to job script..")
            self.jobscript.addMany(_global.host.loaders)
        
        return self.jobscript # return the jobscript to user for debugging
            
    def mdrun(self, new, **kwargs):
        """Execute gmx mdrun"""
        
        def _do(new, **kwargs):
            """Runs mdrun, depending on if jobscript is there or not"""
            if self.jobscript is None:
                self.gmx('mdrun', **kwargs, deffnm = new)
            else:
                self.jobscript.add(self.gmx('mdrun', **kwargs, deffnm = new, onlycmd=True))
        
        if 'noappend' in kwargs: # set noappend if noappend=True passed
            if kwargs['noappend'] == True:
                del kwargs['noappend']
                _do(new, **kwargs, noappend = '')
            else: _do(new, **kwargs)
        else:
            _do(new, **kwargs)
        
        if self.jobscript:
            if 'dirc' in kwargs: # for jobscripts, dirc has to be passed to sbatch() not gmx()
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
            if not _global.host.isLocal():
                _global.logger.write('info', "Please do not close remote host until the job is done!")
    
    def restart(self, new, until=0, extend=0, fromcrash=False, noappend=True):
        """
        Restart a simulation, either:
            - from crash, in which case mdrun is called
            - increase the timesteps of a completed run by extending run time until time t in ps,
              in this case convert-tpr has to called with -until
            - increase timesteps of a completed run by adding time t in ps to current time,
            in this case convert-tpr has to called with -extend
        """
        
        self.setcurrent(new)
        
        if fromcrash:
            out = self.mdrun(new, s = self.getcurrent('tpr'), cpi = self.getcurrent('cpt'), noappend = noappend, dirc=new)
        elif until != 0:
            self.gmx('convert-tpr', s = self.getcurrent('tpr'), until = until, o = f'{new}.tpr', dirc=new)
            out = self.mdrun(new, s = f'{new}.tpr', cpi = self.getcurrent('cpt'), noappend = noappend, dirc=new)
        elif extend != 0:
            self.gmx('convert-tpr', s = self.getcurrent('tpr'), extend = extend, o = f'{new}.tpr', dirc=new)
            out = self.mdrun(new, s = f'{new}.tpr', cpi = '{cpt}.cpt', noappend = noappend, dirc=new)
        
        self.saveToYaml()
        return out
        
        
    def grompp(self, mdp, new, **kwargs):
        """
        Writes mdp to MDP file
        Finds latest gro/trr, and top
        Execute gmx grompp
        """
        new = new.lower().replace(' ', '_')
        # write mdp file into folder dirc
        if 'dirc' in kwargs:
            mdp_file = f'{kwargs["dirc"]}/{new}.mdp'
            # don't delete kwargs['dirc'], it should be passed to gmx
        else:
            mdp_file = f'{new}.mdp'
        _global.host.write(str(mdp), mdp_file)
        
        # a custom gro/trr file can be passed to grompp
        if 'gro' in kwargs:
            gro_file = kwargs['gro']
            del kwargs['gro']
        elif 'trr' in kwargs:
            trr_file = kwargs['trr']
            del kwargs['trr']
        
        # if no custom file is passed, search for it in folders
        gro = self.getcurrentNone(gro_file, 'gro', level=True, exp=False)
        trr = self.getcurrentNone(trr_file, 'trr', level=True, exp=False)
        
        if gro == None and trr ==  None:
            raise MiMiCPyError("No coordinate data (gro/trr file) was found..")
        
        if trr != None: # if both gro and trr were found
            if gro[0] < trr[0]: # find the latest one, and run grompp accordingly
                self.gmx('grompp', f = f'{new}.mdp', c = gro[1],\
                     p = self.getcurrent('top'), o = f"{new}.tpr", **kwargs)
            else:
                self.gmx('grompp', f = f'{new}.mdp', c = gro[1],\
                     p = self.getcurrent('top'), t=trr[1], o = f"{new}.tpr", **kwargs)
        else: # otherwise, just use the gro file
            self.gmx('grompp', f = f'{new}.mdp', c = gro[1],\
                     p = self.getcurrent('top'), o = f"{new}.tpr", **kwargs)
        
    def run(self, mdp, **kwargs):
        """
        Wrapper around gromp() and mdrun()
        Also takes care of directory creation
        """
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
    """
    Runs MiMiC runs by running both gmx mdrun and cpmd
    Inherits from .core.base.BaseHandle
    
    """
    
    # functions to set individual slurm parameters for gmx and cpmd
    def setGMXSettings(self, **kwargs): self.gmx_opt = kwargs
    def setCPMDSettings(self, **kwargs): self.cpmd_opt = kwargs
    
    def run(self, inp, tpr=None, dirc=''):
        """
        Writes inp to CPMD file
        And runs gmx mdrun and cpmd depending
        Also supports Slurm
        """
        new = inp.info.lower().replace(' ', '_')
        self.setcurrent(new)
        _global.host.mkdir(f"{new}/cpmd") # cpmd
        _global.host.mkdir(f"{new}/gmx") # and gmx directories
        
        inp.mimic.paths = f"1\n{_global.host.pwd()+new}/gmx" # set path in cpmd script
        _global.logger.write('debug2', "Set PATH in MIMIC section as {inp.mimic.paths}")
        _global.host.write(str(inp), f"{new}/cpmd/{new}.inp")
        
        tpr = self.getcurrentNone('mimic-tpr') # look for mimic tpr in prepQM folder
        
        _global.host.cp(tpr, f'{new}/gmx/mimic.tpr')
        
        # below is similar procedure to simulate.MD.run()
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
            
            if not _global.host.isLocal():
                _global.logger.write('info', "Please do not close remote host until the job is done!")
            
        self.saveToYaml()
        
        