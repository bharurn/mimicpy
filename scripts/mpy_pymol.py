import mimicpy
from pymol import cmd

# TO DO: Error handling
mimicpy.setLogger(1)

class Main:
    
    @staticmethod
    def load(mpt, gro):
    	Main.qm = mimicpy.Prepare(mpt, gro, selector=mimicpy.selector.PyMOL())
    	
    @staticmethod
    def prepqm(selection, cpmd='cpmd.inp', ndx='index.ndx'):
        Main.qm.add(selection)
        Main.qm.getInp(ndx, cpmd)

cmd.extend('mpyload', Main.load)
cmd.extend('prepqm', Main.prepqm)