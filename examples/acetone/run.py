# example of running mimicpy as a python module

from mimicpy import DefaultSelector, Preparation

# Add topology and structure to the Default selector
selector = DefaultSelector('acetone.top', 'acetone.gro')

# Start preparation session
preparation = Preparation(selector)

# Add atoms to QM partition
preparation.add('resname is ACT')

# Get GROMACS and CPMD inputs
ndx, inp = preparation.get_mimic_input()