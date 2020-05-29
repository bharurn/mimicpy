# MiMiCPy
MiMiCPy is a python package for efficient set-up and execution of QM/MM simulations using the MiMiC CPMD/Gromacs interface developed at Forschungszentrum Juelich and EPFL. [1] It supports a custom selection language and integration with molecular visualization software packages like PyMOL for easy preparation of the MiMiC input scripts. Preparation of Gromacs topologies for non-standard ligands, and other system preparation steps are also supported for the Amber force field. Support for wider system topologies and force fields are in development.

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
The package has been tested and confirmed to work on Mac OSX and Linux. Running on Windows is not supported natively, but can work by using a UNIX shell utility for Windows. However, this has not been tested and performance is unknown.

Currently, at least Python 3.6 is required to run the package. If it is unavailable on the remote server you are using to run MiMiC (and you do not have sudo privileges), the package and required python verions need to only be installed locally, and the MiMiCPy host can be set-up to run all jobs remotely. For more details, please read the docs.

## Demo
Below is a demo for preparing a given protein+ligand system and running a MiMiC simulation:
```python
import mimicpy
prt = mimicpy.loadRCSB('3INM') # load protein, ligands from PDBID 3INM and generate ligand topologies

prep = mimicpy.prepare.MM() # MM prepare handle
prep.getTopology(prt) # generate gromacs, mimicpy topology and solvate

md = mimicpy.simulate.MD.continueFrom(prep) # MD simulation handle
em_mdp = mimicpy.scripts.MDP.defaultEM() # get the default energy minimization MDP Gromacs file
md.run(em_mdp) # minimize system

qm = mimicpy.prepare.QM() # QM prepare handle
qm.add('resname is ICT and id < 13085') # add ligand ICT to the QM region
cpmd_inp = qm.getInp() # get the CPMD input file for a MiMiC run
cpmd_inp.cpmd.molecular__dynamics__cp = '' # set the Car-Parrinello option ON in the input script

mimic = mimicpy.simulate.MiMiC.continuefrom(qm) # MiMiC simulation handle
mimic.run(cpmd_inp) # run a MiMiC simulation using the above cpmd script
```
For more details and options please refer to the docs before implementing MiMiCPy for your research.
 
## References
[1] J. Chem. Theory Comput. 2019, 15, 6, 3810–3823