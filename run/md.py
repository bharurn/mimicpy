from pygmx.run import mdp, gmxrun
from pyshell.remote import slurm

class MD(gmxrun.GMXRun):
    
    def __init__(self, server, job_settings=None, modules=[], sources=[], gmx="gmx", directory=""):
        print("**Setting up MD Engine**")
        
        
        if job_settings:    
            if not job_settings.noCommands():
                print("Warning! Jobscript contains commands..")
            self.slurm_settings = job_settings
        else:
            self.slurm_settings = slurm.Jobscript()
            
        super().__init__(server, modules, sources, gmx, directory)
        
    def _getjob(self, jobscript, name):
        if jobscript:
            return jobscript
        elif self.slurm_settings:
            j = self.slurm_settings
            j.name = name
            return j
        else:
            return slurm.Jobscript()
            
    def restartMD(self, old, new, until=0, extend=0, fromcrash=False, job_settings=None):
        job = self._getjob(job_settings, 'EQMD')
        
        job.source(self.sources)
        job.modules(self.modules)
        
        return self.ssh.restartMD(old, new, job, until=0, extend=0, fromcrash=False)
        
    def _calc(self, mdp_data, old, new, topol, job):
        
        with self.shell.vi(f"{new}.mdp", 'w') as f:
            f.write(str(mdp_data))
        
        job.add(f"{self.gmx} grompp -f {new}.mdp -c {old}.gro -r {old}.gro -p {topol}.top -o {new}.tpr")
        job.add(f"$MPIEXEC $FLAGS_MPI_BATCH {self.gmx} mdrun -deffnm {new}")
        
        return job
    
    def _run_job(self, job):
        
        jid = self.shell.sbatch(job)
        
        print(f"Job started with job ID: {jid}..\nPlease use squeue() to check status..\n**Done**")
        
        return jid
    
    def calc(self, mdp_data, old, new, topol='topol', job_settings=None):
        
        print("**Running MD Simulation**")
        print("Setting up Slurm job..")
        
        job = self._getjob(job_settings, 'EM')
            
        job.source(self.sources)
        job.module(self.modules)
        
        job = self._calc(mdp_data, old, new, topol, job_settings)
        
        return self._run_job(job)
    
    def _em(self, job):
        
        job.add('$MPIEXEC $FLAGS_MPI_BATCH gmx_mpi mdrun -deffnm em')
        
        return job
    
    def em(self, job_settings=None):
        
        print("**Running Energy Minimization**")
        print("Setting up Slurm job..")
        
        job = self._getjob(job_settings, 'EM')
            
        job.source(self.sources)
        job.module(self.modules)
        
        job = self._em(job)
        
        return self._run_job(job)
        
    def nvt(self, nvt_mdp=mdp.MDP.defaultNVT(), old='em', new='nvt', topol='topol', job_settings=None):
        
        print("**Running NVT Simulation**")
        print("Setting up Slurm job..")
        
        job = self._getjob(job_settings, new.upper())
            
        job.source(self.sources)
        job.module(self.modules)
        
        job = self._calc(nvt_mdp, old, new, topol, job)
        
        return self._run_job(job)
    
    def npt(self, npt_mdp=mdp.MDP.defaultNPT(), old='nvt', new='npt', topol='topol', job_settings=None):
        
        print("**Running NPT Simulation**")
        print("Setting up Slurm job..")
        
        job = self._getjob(job_settings, new.upper())
            
        job.source(self.sources)
        job.module(self.modules)
        
        job = self._calc(npt_mdp, old, new, topol, job)
        
        return self._run_job(job)
    
    def eqmd(self, eq_mdp=mdp.MDP.defaultEQMD(), old='npt', new='md', topol='topol', job_settings=None):
        
        print("**Running Production MD**")
        print("Setting up Slurm job..")
        
        job = self._getjob(job_settings, new.upper())
            
        job.source(self.sources)
        job.module(self.modules)
        
        job = self._calc(eq_mdp, old, new, topol, job)
        
        return self._run_job(job)
            