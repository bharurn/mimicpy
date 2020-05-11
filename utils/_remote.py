import time
#import re

class Shell:
    def __init__(self, name, config):
        
        self.name = name
        
        self.config = config
        
        conf = self._lookup(name)
        
        if 'proxyjump' in conf:
            proxy_name, port = conf['proxyjump'].split(':')
            vmtransport = self._getssh(proxy_name).get_transport()
            dest_addr = (self._lookup(name)['hostname'], int(port))
            local_addr = (self._lookup(proxy_name)['hostname'], int(port))
            vmchannel = vmtransport.open_channel("direct-tcpip", dest_addr, local_addr)
            
            self.ssh = self._getssh(name, sock=vmchannel)
        
        else:
            self.ssh =  self._getssh(name)
            
        self.sftp = self.ssh.open_sftp()
        self.shell = self.ssh.invoke_shell()
        self.stdout = ''
        self.dir = '.'
        self.parallel = False
        self.errorHandle = None
        self.stdout_source = '_std'
        self.nohup = False
        self.decoder = 'utf-8'
    
    def _lookup(self, name):
        try:
            import paramiko
        except ImportError:
            raise Exception("Install paramiko python package to remotely run MiMiC")
            
        config = paramiko.SSHConfig()
        config.parse(open(self.config))
        conf = config.lookup(name)
        
        return conf
    
    def _getssh(self, name, sock=None):
        try:
            import paramiko
        except ImportError:
            raise Exception("Install paramiko python package to remotely run MiMiC")
            
        conf = self._lookup(name)
    
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
    
    def run(self, cmd, stdin=None, errorHandle=None, onNewChan=None, query_rate=0.3):
        
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
                
            #if 'mdrun' in cmd: log = True
            #else: log = False
        
            prev = startout
            while not self.parallel:
                if query_rate > 0:    
                    time.sleep(query_rate)
                
                self.queryStdout()
                
                if prev != self.stdout:
                    text = self.stdout.replace(prev, '')
                    
                    prev = self.stdout
                    
                    #if log:
                        #x = re.compile(r"^\s*(Step\s*Time)\n(.+)",\
                        #               re.MULTILINE)
                        #res = x.findall(text)
                        #if res:
                        #    for r in res: print(f'Step: {r[1].split()[0]}   |   Time: {r[1].split()[1]}')
                    
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
    
    def clearStdout(self):
        self.stdout = ''
    
    def __enter__(self):
        return self
        
    def __exit__(self, exc_type, exc_val, exc_tb):
        self.ssh.close()