from pygmx.run import ssh

class GMXRun:
    def __init__(self, server, modules=[], sources=[], gmx="gmx", directory=".", ignloadxerr = True):
        print(f"Using {gmx} as Gromacs executable..")
        self.gmx = gmx
        print(f"Invoking interactive SSH connection to {server}..")
        self.shell = ssh.GmxSSH(server, self.gmx, directory)
        
        self.modules = modules
        
        for module in modules:
            print(f"Loading module {module}..")
            try:
                self.shell.run(f'module load {module}')
            except Exception as e:
                if ignloadxerr:
                    print(f"Warning! {e}. Ignoring..")
                else:
                    raise Exception(str(e))
        
        self.sources = sources
        
        for source in sources:
            print(f"Sourcing {source}..")
            try:
                self.shell.run(f'source {source}')
            except Exception as e:
                if ignloadxerr:
                    print(f"Warning! {e}. Ignoring..")
                else:
                    raise Exception(str(e))
        
        print('**Done**')
        
    @classmethod
    def continueSession(cls, session):
        return cls(session.server, session.sources, session.gmx, session.directory)
    
    def saveSessionToFile(self, filename):
        pass
    
    @classmethod
    def continueSessionFromFile(cls, filename):
        pass
    
    def restartShell(self, new_shell, fullrestore=True):
        self.shell.close()
        self.shell = new_shell
        
        if fullrestore:
            self.shell.cd(self.pwd(), mkdir=True)
        
            for module in self.modules:
                self.shell.run(f'module load {module}')
        
            for source in self.sources:
                self.shell.run(f'source {source}')
        else:
            self.shell.dir = ''
            self.modules = ''
            self.sources = ''
    
    def __enter__(self):
        return self
        
    def __exit__(self, exc_type, exc_val, exc_tb):
        self.__del__()
    
    def __del__(self):
        self.shell.__del__()