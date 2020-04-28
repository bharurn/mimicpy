from prettytable import PrettyTable

class Jobscript:
    def __init__(self, name='JOB', shebang='/bin/bash', allocation='', cmds = [], **kwargs):
        self._shebang = shebang
        self._allocation = allocation
        self._name = name
        self._dir = ""
        self.nodes = 2
        self.ntasks = 16
        self.cpus_per_task = 3
        self.mem_per_cpu = '1240M'
        self.export = 'NONE'
        self.time = '00:10:00'
        
        for key, value in kwargs.items():
            setattr(self, key, value)
        
        self._cmds = cmds
    
    def setDir(self, dirc):
        self._dir = f"./{dirc}"
        
    @classmethod    
    def loadFromPrevious(cls, script):
        if isinstance(script, str):
            text = script
        else:
            text = script.read().decode('utf-8')
            
        shebang = "/usr/sh"
        name = "JOB"
        cmds = []
        kwargs = {}
        for line in text.splitlines():
            if line.startswith('#!'):
                  shebang = line[2:].strip()
            elif line.startswith('#SBATCH -A'):
                  allocation = line[11:].strip()  
            elif line.startswith('#SBATCH --'):
                l = line[10:].replace('-', '_').split('=')
                key = l[0].strip()
                val = l[1].strip()
                if key == 'job_name':
                    name = val
                elif key == 'output':
                    continue
                elif key != 'name':
                    kwargs.update({key:val})
            elif line.strip() != '':
                cmds.append(line)
        
        return cls(name=name, shebang=shebang, allocation=allocation, cmds=cmds, **kwargs)
                
    def add(self, cmd):
        self._cmds.append(cmd)
    
    def addMany(self, cmds):
        self._cmds.extend(cmds)
    
    def source(self, sources):
        for source in sources:
            self._cmds.append(f'source {source}')
        
    def module(self, modules):
        for module in modules:
            self._cmds.append(f'module load {module}')
    
    def __getattr__(self, val):
        if val == 'job_name' or val == 'name':
            return self._name
        elif val == 'output':
            return f'{self._name}.%J.out'
    
    def noCommands(self):
        if len(self._cmds) <= 0:
            return True
        return False
    
    def clearCommands(self, keep_source=False, keep_module=False):
        if keep_source:
            s = []
            for i, c in enumerate(self._cmds):
                if c.split()[0] == 'source':
                    s.append(self._cmds[i])
            self._cmds = s.copy() 
        elif keep_module:
            s = []
            for i, c in enumerate(self._cmds):
                if c.split()[0] == 'module':
                    s.append(self._cmds[i])
            self._cmds = s.copy()
        else:            
            self._cmds = []
        
    def __str__(self):
        cmd = f'#!{self._shebang}\n'
        
        data = [attr for attr in dir(self) if not callable(getattr(self, attr)) and not attr.startswith('_')]
        
        for d in data:
            if getattr(self,d) == None: continue
            d_ = d.replace('_','-')
            cmd += f"#SBATCH --{d_}={getattr(self, d)}\n"
        
        cmd += f'#SBATCH --job-name={self._name}\n'
        cmd += f'#SBATCH --output={self._name}.%J.out\n'
        if self._allocation:
            cmd += f'#SBATCH -A {self._allocation}\n\n'
        
        cmd += '\n'.join(self._cmds)
        
        return cmd
    
class Jobinfo:
    def __init__(self, jobid, partition, name, user, status, time, nodes, nodelist):
        self.jobid = jobid
        self.partition = partition
        self.name = name
        self.user = user
        self.status = status
        self.time = time
        self.nodes = nodes
        self.nodelist = nodelist
        
    def __str__(self):
        return str(self.tolist())
    
    def __len__(self):
        return 8
    
    def __iter__(self, i):
        if i == 0: return self.jobid
        elif i==1: return self.partition
        elif i==2: return self.name
        elif i==3: return self.user
        elif i==4: return self.status
        elif i==5: return self.time
        elif i==6: return self.nodes
        elif i==7: return self.nodelist
    
    def tolist(self):
        s = []
        for i in range(8):
            s.append(self.__iter__(i))
        
        return s
        
class Joblist:
    def __init__(self, jobs, ids, names):
        self.jobs = jobs
        self.jobid = dict(zip(ids, self.jobs))
        self.jobname = dict(zip(names, self.jobs))
    
    @classmethod
    def fromQueue(cls, q):
        jobs = []
        jobid = []
        jobname = []
        
        if len(q.splitlines()) == 1:
            raise Exception("Error: No jobs in queue!")
            
        for i in q.splitlines()[1:]:
            a = i.split()
            jobs.append(Jobinfo(*a))
            jobid.append(int(a[0]))
            jobname.append(a[2])
        
        return cls(jobs, jobid, jobname)
    
    def __str__(self):
        t = PrettyTable(["JobID", "Partition", "Name", "User", "Status", "Time", "Nodes", "Nodelist"])
        for i in self.jobs:
            t.add_row(i.tolist())
        
        return t.get_string()
    
    def __getitem__(self, i):
        if isinstance(i, int):
            return self.jobid[i]
        else:
            return self.jobname[i]
