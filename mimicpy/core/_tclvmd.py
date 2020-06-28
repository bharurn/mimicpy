from ..utils.errors import MiMiCPyError

class CustomAtomSel:
    """Class to mimic AtomSel class from vmd python module"""
    
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

class TclVMDConnector:
    """
    Class to mimic vmd python module
    Results of following Tcl commands to be passed in the constructor:
        ##Values set to TclVMDConnector.sele
        atomsel get name
        atomsel get type
        atomsel get index
        atomsel get mass
        atomsel get element
        atomsel get resname
        atomsel get resid
        atomsel get x
        atomsel get y
        atomsel get z
        ##
        ##Values set to TclVMDConnector.box_size
        molinfo <molid> get a
        molinfo <molid> get b
        molinfo <molid> get c
        molinfo <molid> get alpha
        molinfo <molid> get beta
        molinfo <molid> get gamma
        ##
    All should be passed in this order, each param being a space seperated strings (this is Tcl's default behaviour)
    The gro file should be loaded in Tcl/VMD, it will not be loaded here
    """
    
    def __init__(self, params):
        if len(params) < 16:
            raise MiMiCPyError("Not enough params passed to TclVMDConnector")
            
        self.sele = params[:10]
        self.box_size = params[10:]
        # set props of molecules, molecule.load() -> return dummy molid, molecule.get_periodc() -> return actual box size
        self.molecule = type('obj', (object,), {'load' : lambda a,b: -1, 'get_periodic': self.__get_periodic})
        
    def __get_periodic(self, a=-1, b=-1):
        if not self.box_size:
            raise MiMiCPyError("Did not receive system box size information from Tcl")
        
        box_vals = [float(b) for b in self.box_size]
        box_keys = ['a', 'b', 'c', 'alpha', 'beta', 'gamma']
        
        # return as dict of vals like vmd python module
        return dict(zip(box_keys, box_vals))
        
        
    def atomsel(self, selection, molid):
        # atomsel method to return CustomAtomSel object
        if not self.sele:
            raise MiMiCPyError("Did not receive QM atoms information from Tcl")
        
        # get selection params from tcl script
        params = ['name', 'type', 'index', 'mass', 'element', 'resname', 'resid', 'x', 'y', 'z']
        # self.sele if expected to be list of strings of name, type, resid,.. directly from tcl script
        kwargs = dict(zip(params, self.sele))
        
        return CustomAtomSel(**kwargs)