# MiMiCPy
MiMiCPy is the python library for preparing and executing QM/MM simulations with the MiMiC CPMD/Gromacs interface developed at Forschungszentrum Juelich and EPFL. [1] For more details on compiling the MiMiC source code, please refer to [2].

MiMiCPy comes with a set of command lines tools to prepare MiMiC input scripts, and integration with molecular visualization packages. These features, as well as running MiMiC simulations, can also be accessed through the python interface.

## Installation
MiMiCPy is not available on pip/conda yet. To install run the following command in the terminal:
```
git clone https://github.com/bharurn/mimicpy
```
and and run ```setup.py``` .

## Portability Issues
The package has been tested and confirmed to work on Linux and MacOS systems. Running on Windows is not fully supported and has not been tested. Currently, at least Python 3.6 is required to run the package.

## Demo
Below is a demo for setting up CPMD and GROMACS input scripts for a MiMiC simulation:
```python
import mimicpy

```
For more details and options please refer to the docs.
 
## References
[1] J. Chem. Theory Comput. 2019, 15, 6, 3810–3823

[2] MiMiC Hybrid Quantum Mechanical Interface (http://manual.gromacs.org/documentation/2019-rc1/reference-manual/special/mimic-qmmm.html)
