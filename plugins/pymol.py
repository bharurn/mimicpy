
##### MiMiCPy PyMOL settings script
##
import mimicpy
from pymol import cmd

def prepqm(mpt, selection=None, is_link=False, cpmd='cpmd.inp', ndx='index.ndx'):
    qm = mimicpy.Preparation(mimicpy.PyMOL(mpt))
    qm.add(selection, is_link)
    qm.get_mimic_input(ndx_out=ndx, inp_out=cpmd)

cmd.extend('prepqm', prepqm)
##
##################################