# MiMiCPy
MiMiCPy is a python package for efficient set-up and execution of QM/MM simulations using the newly developed MiMiC QM/MM interface developed at Forschungszentrum Juelich and EPFL. The package is currently optimized for simulating protein (including handling non standard ligands). It can also run pure MM simulation using Gromacs back-end, and has support for submitting jobs to a remote server. It is currently in alpha testing, and more features will be added soon.

## Installation
MiMiCPy is not available on pip/conda yet. To install run the following command in the terminal:
```
git clone https://github.com/bharurn/mimicpy
```
and include this as the first line of your Python script:
```python
import sys; sys.path.insert(0,'/home/usr/mimicpy/')
```
where `/home/usr/mimicpy/` is the path to the cloned repository. The PYTHONPATH enviornment variable can also be edited to permanently affect this change.

For the required packages/dependencies please read the requirements.txt and manually install all necessary packages.
These issues are only a temporary hassle, and will be fixed once the package is stable enough to publish it on pip/conda.

## Portability Issues
The package has only been tested on Mac OSX and Linux, and it is not expected to be stable on Windows machines. Currently, at least Python 3.6 is required to run the package. If it is unavailable on the remote server you are using to run MiMiC (and you do not have sudo privileges), the package and required python verions need to only be installed locally, and the MiMiCPy host can be set-up to run all jobs remotely. For more details, please read the docs.

## Demo
Below is a demo for preparing a given protein+ligand system and running a MiMiC simulation:
```python
import mimicpy
prt = mimicpy.loadRCSB('3INM') # load protein, ligands from PDBID and generate topologies

prep = mimicpy.prepare.MM() # MM prepare handle
prep.getTopology(prt) # generate gromacs topology of protein+ligand system
prep.getBox() # generate solavted box and neutralize charge with ions

md = mimicpy.simulate.MD.continueFrom(prep) # MD simulation handle
em_mdp = mimicpy.scripts.MDP.defaultEM() # get the default energy minimization MDP Gromacs file
md.run(em_mdp) # minimize system

qm = mimicpy.prepare.QM() # QM prepare handle
qm.add('resName ICT and chainID A') # add ligand ICT chain A to the QM region
cpmd_inp = qm.getInp() # get the CPMD input file for a MiMiC run
cpmd_inp.cpmd.molecular__dynamics__cp = '' # set the Car-Parrinello option ON in the input script

mimic = mimicpy.simulate.MiMiC.continuefrom(qm) # MiMiC simulation handle
mimic.run(cpmd_inp) # run a MiMiC simulation using the above cpmd script
```
For more details and options please refer to the docs.
