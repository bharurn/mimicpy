import time
from stat import S_ISDIR, S_ISREG

class Shell:
    def __init__(self, name, config):
        
        self.name = name
        
        self.config = config
        
        conf = self._lookup(name)
        
        if 'proxyjump' in conf:
            proxy_name, port = conf['proxyjump'].split(':')
            vmtransport = Shell._getssh(proxy_name).get_transport()
            dest_addr = (Shell._lookup(name)['hostname'], int(port))
            local_addr = (Shell._lookup(proxy_name)['hostname'], int(port))
            vmchannel = vmtransport.open_channel("direct-tcpip", dest_addr, local_addr)
            
            self.ssh = Shell._getssh(name, sock=vmchannel)
        
        else:
            self.ssh =  Shell._getssh(name)
            
        self.sftp = self.ssh.open_sftp()
        self.shell = self.ssh.invoke_shell()
        self.stdout = ''
        self.dir = '.'
        self.query_rate = 0.3
        self.parallel = False
        self.errorHandle = None
        self.stdout_source = '_std'
        self.nohup = False
        self.decoder = 'utf-8'
    
    def module(self, module):
        self.run(f'module load {module}')
        
    def source(self, source):
        self.run(f'source {source}')
        
    @staticmethod
    def _lookup(self, name):
        try:
            import paramiko
        except ImportError:
            raise Exception("Install paramiko python package to remotely run MiMiC")
            
        config = paramiko.SSHConfig()
        config.parse(open(self.config))
        conf = config.lookup(name)
        
        return conf
    
    @staticmethod
    def _getssh(name, sock=None):
        try:
            import paramiko
        except ImportError:
            raise Exception("Install paramiko python package to remotely run MiMiC")
            
        conf = Shell._lookup(name)
    
        ssh = paramiko.SSHClient()
        ssh.set_missing_host_key_policy(paramiko.AutoAddPolicy())
        if sock is None:
            ssh.connect(hostname=conf['hostname'], username=conf['user'])
        else:
            ssh.connect(hostname=conf['hostname'], username=conf['user'], sock=sock)
        
        return ssh
    
    def __del__(self): self.ssh.close()
        
    def queryStdout(self, buff = 1024):
        
        if self.stdout_source != '_std':
            try:
                self.sftp.stat(self.stdout_source)
            except:
                self.stdout = ''
                return
            
            with self.sftp.open(self.stdout_source, 'r') as f:
                self.stdout = f.read().decode(self.decoder)
        else:        
            while self.shell.recv_ready():
                self.stdout += self.shell.recv(buff).decode(self.decoder)
        
    def redirectStdout(self, fname):
        self.stdout_source = fname
    
    def run(self, cmd, stdin=None, errorHandle=None, onNewChan=None):
        
        if errorHandle is None: errorHandle = self.errorHandle
        
        if onNewChan:
            sin, out, err = self.ssh.exec_command(cmd)
            
            if err.read().decode(self.decoder).strip() != '':
                if errorHandle:
                    errorHandle(err.read().decode(self.decoder).strip())
                else:
                    raise Exception(f'Error! {err.read().decode(self.decoder)} {out.read().decode(self.decoder)}')
                
            return out.read().decode(self.decoder)
        
        else:
            self.queryStdout()
        
            startout = self.stdout
            
            if stdin is None:
                self.shell.send(cmd+'\n')
            else:
                self.shell.send(f'printf "{stdin}" | {cmd}\n')
        
            prev = ''
            while not self.parallel:
                if self.query_rate > 0:    
                    time.sleep(self.query_rate)
                
                self.queryStdout()
                
                if prev != self.stdout:
                    prev = self.stdout
                    
                    text = self.stdout.replace(startout, '')
                    
                    if errorHandle != None:
                        errorHandle(text)
                    else:
                        if 'ERROR' in text.upper():
                            for line in self.stdout.splitlines()[::-1]:
                                if 'ERROR' in line.upper():
                                    raise Exception(line)
                else:
                    break
        
            lines = self.stdout.replace(startout, '').splitlines()[1:-1]
        
            return '\n'.join(lines)
    
    def sbatch(self, job):
        if job.noCommands():
            raise Exception("No commands in jobscript!")
        job.setDir(self.dir)
        
        f = self.vi(f'{job.name}.sh', 'w')
        f.write(str(job))
        f.close()
        
        def _sbatch_err(err):
            if 'error' in err.lower():
                raise Exception(err) 
            
        out = self.run(f'cd {self.dir} && sbatch {job.name}.sh', errorHandle=_sbatch_err, onNewChan=True)
        
        try:
            idx = int(out.split()[3])
        except:
            raise Exception(out)
        
        return idx
        
    def cd(self, directory, mkdir=False):
        if directory.strip() == '.':
            return -1
        
        if not self.fileExists(directory):
             if mkdir:       
                self.sftp.mkdir(directory)
                return 1
             else:
                 raise Exception(f'Directory {directory} not found')
        else:     
            self.dir = directory+'/'
            self.run(f'cd {self.dir}')
            self.sftp.chdir(self.dir)
            return 0
    
    def scancel(self, jobid):
        return self.run(f'scancel {jobid}', onNewChan=True)
    
    def vi(self, name, s):
        return self.sftp.open(f"{name}", s)
    
    def clearStdout(self):
        self.stdout = ''
    
    def fileExists(self, file):
        try:
            self.sftp.stat(file)
            return True
        except FileNotFoundError:
           return False
    
    def ls(self, file_eval=lambda a: True, dir_eval=lambda a: True):
        
        files = []
        
        for entry in self.sftp.listdir_attr():
            if S_ISDIR(entry.st_mode) and dir_eval(entry.filename):
                files.append(entry.filename)
            elif S_ISREG(entry.st_mode) and file_eval(entry.filename):
                files.append(entry.filename)
            
        return files
    
    def pwd(self): return self.dir
    
    def put(self, local, remote): self.sftp.put(local, remote, confirm=True) #to do use callback to create a tqdm loader
    
    def get(self, remote, local): self.sftp.get(remote, local) #to do use callback to create a tqdm loader
    
    def __enter__(self):
        return self
        
    def __exit__(self, exc_type, exc_val, exc_tb):
        self.ssh.close()