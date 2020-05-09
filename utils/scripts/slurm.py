from .base import Script

class Slurm(Script):
    def __init__(self, name='JOB', shebang='/bin/bash', allocation='', cmds = [], **kwargs):
        if kwargs is None:
            kwargs = {'nodes':2, 'ntasks': 16, 'cpus_per_task': 3, 'mem_per_cpus': 3, 'export': 'NONE', 'time': '00:10:00'}
            
        super().__init__(**kwargs)
        
        self._shebang = shebang
        self._allocation = allocation
        self._name = name
        self._dir = ""
        self._cmds = cmds
    
    def setDir(self, dirc):
        self._dir = f"./{dirc}"
    
    def setSpecial(self, name=None, shebang=None, allocation=None):
        if name: self._name = name
        if shebang: self._shebang = shebang
        if allocation: self._allocation = allocation
    
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
        else:
            return super().__getattr__(val)
    
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
        
        for d in self.params():
            if getattr(self,d) == None: continue
            d_ = d.replace('_','-')
            cmd += f"#SBATCH --{d_}={getattr(self, d)}\n"
        
        cmd += f'#SBATCH --job-name={self._name}\n'
        cmd += f'#SBATCH --output={self._name}.%J.out\n'
        
        if self._allocation:
            cmd += f'#SBATCH -A {self._allocation}\n\n'
        
        cmd += '\n'.join(self._cmds)
        
        return cmd