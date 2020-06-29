# MiMiCPy
MiMiCPy is the python library for preparing and executing QM/MM simulations with the MiMiC CPMD/Gromacs interface developed at Forschungszentrum Juelich and EPFL. [1] For more details on compiling the MiMiC source code, please refer to [2].

MiMiCPy comes with a set of command lines tools to prepare MiMiC input scripts, and integration with molecular visualization packages. These features, as well as running MiMiC simulations, can also be accessed through the python interface.

## Installation
MiMiCPy is not available on pip/conda yet. To install run the following command in the terminal:
```
git clone https://github.com/bharurn/mimicpy
```
and include this as the first line of your Python script:
```python
import sys; sys.path.insert(0,'/path/to/mimicpy/')
```
where `/path/to/mimicpy/` is the path to the cloned repository. The PYTHONPATH environment variable can also be edited to permanently affect this change.

For the required packages/dependencies please read the requirements.txt and manually install all necessary packages.
These difficulties in setup are only a temporary hassle, and will be fixed once the package is stable enough to be published on pip/conda.

## Portability Issues
The package has been tested and confirmed to work on Linux and MacOS systems. Running on Windows is not fully supported and has not been tested. Currently, at least Python 3.6 is required on the local machine used to run the package.

## Demo
Below is a demo for running a MiMiC simulation:
```python
import mimicpy

md = mimicpy.simulate.MD('topol_dir') # MD simulation handle
em_mdp = mimicpy.scripts.MDP.defaultEM() # get the default energy minimization MDP Gromacs file
md.run(em_mdp) # minimize system

mpt = mimicpy.parsers.MPT.fromTop(topol) # get MPT file
qm = mimicpy.Prepare(mpt, 'conf1.gro') # QM prepare handle, read MPT file
qm.add('resname is GLY and resid 450') # add ligand GLY to the QM region
cpmd_inp = qm.getInp() # get the CPMD input file for a MiMiC run

cpmd_inp.cpmd.molecular__dynamics__cp = '' # set the Car-Parrinello option ON in the input script
mimic_mdp = mimicpy.scripts.MDP.defaultMiMiC() # get the default MiMiC MDP Gromacs file
mimic = mimicpy.simulate.MiMiC() # MiMiC simulation handle
mimic.run(mimic_mdp, cpmd_inp) # run a MiMiC simulation using the above cpmd script
```
For more details and options please refer to the docs.
 
## References
[1] J. Chem. Theory Comput. 2019, 15, 6, 3810–3823

[2] MiMiC Hybrid Quantum Mechanical Interface (http://manual.gromacs.org/documentation/2019-rc1/reference-manual/special/mimic-qmmm.html)
