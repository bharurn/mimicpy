class Base():
    
    def __init__(self, modules=[], sources=[], directory=".", loaders=[], ignloaderr = True):
        
        if directory.strip() == '':
            directory = '.'
        
        ret = self.cd(directory, mkdir=True)
        
        if ret == -1:
            pass
        else:
            if ret == 1:
                print(f"{directory} not found, creating new directory..")
            
            print(f"Setting current directory to {directory}..")
            
        self.modules = modules
         
        for module in modules:
            print(f"Loading module {module}..")
            try:
                self.run(f'module load {module}')
            except Exception as e:
                if ignloaderr:
                    print(f"Warning! {e}. Ignoring..")
                else:
                    raise Exception(str(e))
        
        self.sources = sources
        
        for source in sources:
            print(f"Sourcing {source}..")
            try:
                self.run(f'source {source}')
            except Exception as e:
                if ignloaderr:
                    print(f"Warning! {e}. Ignoring..")
                else:
                    raise Exception(str(e))
        
        self.loaders = loaders
        
        for loader in loaders:
            print(f"Running commands: {loader}..")
            try:
                self.run(loader)
            except Exception as e:
                if ignloaderr:
                    print(f"Warning! {e}. Ignoring..")
                else:
                    raise Exception(str(e))
      
    def source(self, sources):
         for source in sources:
            self.run(f'source {source}')
        
    def module(self, modules):
        for module in modules:
            self.run(f'module load {module}')