from pygmx.run import gmxrun

class Analysis(gmxrun.GMXRun):
    
    def _parse(self, out):
        
        with self.shell.vi(out, 'r') as f:
            output = f.read().decode('utf-8')
        
        return self.shell._notes(output), self.shell._errorHandle(output, dont_raise=True)
    
    def parseOutput(self, file_eval=None):
        
        if file_eval is None:
            files = self.shell.ls(file_eval=lambda a: True if a.endswith('.log') or a.endswith('.out') else False)
        else:
            files = self.shell.ls(file_eval=file_eval)
        
        for file in files:
            print(f"\n<========{file}========>")
            
            notes, errors = self._parse(file)
            
            if notes: print(f"\nThere were note(s)/warning(s):\n\n{notes}")
            if errors: print(f"\nThere were error(s):\n\n{errors}")