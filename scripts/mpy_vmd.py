import sys
import pandas as pd
import mimicpy

#TO DO: error handling

mimicpy.setLogger(1)

class CustomAtomSel:
    """Class to mimic AtomSel class from vmd-python"""
    
    def __init__(self, **kwargs):
        # assign name, type, resid, resname, etc. to this class
        # this will be accessed in selector.VMD()
        for k, v in kwargs.items():
            # tcl returns eveything as space separated strings
            lst = v.split()
            
            # convert string to int/float
            if k == 'index' or k == 'resid':
                lst = [int(i) for i in lst]
            elif k == 'mass' or k == 'x' or k == 'y' or k == 'z':
                lst = [float(i) for i in lst]
            setattr(self, k, lst)

class VMDHandle:
    """Class to mimic vmd module"""
    
    def __init__(self):
        self.sele = None
        
    def atomsel(self, selection, molid):
        # atomsel method to return CustomAtomSel object
        if self.sele is None:
            print("Please select the QM atoms")
            sys.exit(1)
        
        return self.sele
        
        

if __name__ == '__main__':
    mpt = sys.argv[1]
    gro = sys.argv[2]
    cpmd = sys.argv[3]
    ndx = sys.argv[4]
    # sys.argv[5:15] is name, type, resid, etc..
       
    qm = mimicpy.Prepare(mpt, gro, selector=mimicpy.selector.VMD(vmd_handle=VMDHandle(), load=False))
        
    # get selection params from tcl script
    params = ['name', 'type', 'index', 'mass', 'element', 'resname', 'resid', 'x', 'y', 'z']
    kwargs = dict(zip(params, sys.argv[5:15]))
        
    qm.selector.cmd.sele = CustomAtomSel(**kwargs) # set the CustomAtomSel object from params
        
    qm.add("sele") # selection doesn't matter
    qm.getInp(cpmd, ndx)