import logging
from ._version import __version__
from ._authors import __authors__
from mimicpy.core.prepare import Preparation
from mimicpy.core.selector import DefaultSelector, VMD, PyMOL
from mimicpy.topology.mpt import Mpt
from mimicpy.topology.top import Top
from mimicpy.scripts.mdp import Mdp
from mimicpy.scripts.cpmd import CpmdScript
from mimicpy.scripts.ndx import Ndx
from mimicpy.coords.base import CoordsIO
from mimicpy.coords.gro import Gro
from mimicpy.coords.pdb import Pdb

logging.basicConfig(# filename='mimicpy.log', format='%(asctime)s %(levelname)s %(message)s',
                    format='%(message)s',
                    filemode='w',
                    level=logging.INFO)
