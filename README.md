# MiMiCPy
MiMiCPy is the python library for preparing and executing QM/MM simulations with the MiMiC CPMD/Gromacs interface developed at Forschungszentrum Juelich and EPFL. [1] For more details on compiling the MiMiC source code, please refer [2].

MiMiCPy is primarily a python wrapper around MiMiC for efficient running of simulations, and preparing input files with support for custom selection language and integration with molecular visualization software packages like PyMOL. It also includes the MiMiCPy CLI tools, a collection of command line tools to the run the core functionalities of the library through the command line.

## Installation
MiMiCPy is not available on pip/conda yet. To install run the following command in the terminal:
```
git clone https://github.com/bharurn/mimicpy
```
and include this as the first line of your Python script:
```python
import sys; sys.path.insert(0,'/path/to/mimicpy/')
```
where `/path/to/mimicpy/` is the path to the cloned repository. The PYTHONPATH enviornment variable can also be edited to permanently affect this change.

For the required packages/dependencies please read the requirements.txt and manually install all necessary packages.
These difficulties in setup are only a temporary hassle, and will be fixed once the package is stable enough to be published on pip/conda.

## Portability Issues
The package has been tested and confirmed to work on Mac OSX and Linux. Running on Windows is not fully supported, however running MiMiCPy on a Windows machine and submitting to a Unix remote host with MiMiC should work. However, this has not been tested.

Currently, at least Python 3.6 is required to run the package. If it is unavailable on the remote server you are using to run MiMiC (and you do not have sudo privileges), MiMiCPy and the dependencies need to only be installed locally, and `mimicpy.setHost()` can be used to set-up all jobs remotely. For more details, please read the docs.

## Demo
Below is a demo for running a MiMiC simulation for a protein-ligand system:
```python
import mimicpy

prep = mimicpy.prepare.MM('topol_data') # feed gromacs topology and coords in topol_data/ to MM prepare handle
prep.addLig('ICT', 'conf1.pdb', 'ligands.itp') # include non-standard ligand ICT
prep.getMPT() # write mimicpy topology

md = mimicpy.simulate.MD.continueFrom(prep) # MD simulation handle
em_mdp = mimicpy.scripts.MDP.defaultEM() # get the default energy minimization MDP Gromacs file
md.run(em_mdp) # minimize system

qm = mimicpy.prepare.QM() # QM prepare handle, read MPT file
qm.add('resname is ICT and resid is 832') # add ligand ICT, chain A to the QM region
cpmd_inp = qm.getInp() # get the CPMD input file for a MiMiC run
cpmd_inp.cpmd.molecular__dynamics__cp = '' # set the Car-Parrinello option ON in the input script

mimic = mimicpy.simulate.MiMiC.continuefrom(qm) # MiMiC simulation handle
mimic.run(cpmd_inp) # run a MiMiC simulation using the above cpmd script
```
For more details and options please refer to the docs.
 
## References
[1] J. Chem. Theory Comput. 2019, 15, 6, 3810–3823

[2] MiMiC Hybrid Quantum Mechanical Interface (http://manual.gromacs.org/documentation/2019-rc1/reference-manual/special/mimic-qmmm.html)
