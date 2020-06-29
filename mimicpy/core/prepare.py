# This is part of MiMiCPy

"""

This module contains the handles to prepare the
MM topology, MPT and QM region/CPMD script

"""

from .selector import Selector
from .._global import _Global as _global
from . import _qmhelper
from ..parsers.mpt import MPT
from ..utils.constants import bohr_rad
from ..scripts import cpmd
from ..utils.errors import MiMiCPyError
import pandas as pd

class Prepare:
    """
    Prepare QM input script, by reading MPT file
    
    """
    
    def __init__(self, mpt, gro, selector=Selector()):
        """Class constructor"""
        
        if isinstance(mpt, MPT): self.mpt = mpt
        else: self.mpt = MPT.fromFile(mpt)
        self.selector = selector
        
        # load mpt and gro into selector, returns box size in gro
        self._mm_box = self.selector.load(self.mpt, gro)
        
        self.qmatoms = pd.DataFrame()
    
    def add(self, selection=None, link=False):
        """Add dataframe to QM region"""
        
        qdf = _qmhelper.cleanqdf( self.selector.select(selection) )
        
        # add a new column link which is integer of link argument
        qdf.insert(2, 'link', [int(link)]*len(qdf))
        
        # add qdf to self.qmatoms, append if already exists
        self.qmatoms = self.qmatoms.append(qdf)
    
    def delete(self, selection=None):
        """Delete from QM region using selection langauage"""
        qdf = _qmhelper.cleanqdf( self.selector.select(selection) )
        # drop selection, ignore pandas errors to disregard extra residues selected that are not present in self.qmatoms
        self.qmatoms = self.qmatoms.drop(qdf.index, errors='ignore')
    
    def clear(self):
        """Empty the QM region to start over"""
        self.qmatoms = pd.DataFrame()
    
    def getInp(self, ndx=None, inp_file=None):
        """
        Create the QM region from the atoms added
        Steps done:
            Assign QM atoms as list to inp._ndx
            If named fiven, writes index files of atoms in QM region with name QMatoms
            Sort self.qmatoms by link and elements
            Read element symbols and all charge info from self.qmatoms
            Add MIMIC, SYSTEM, DFT sections of CPMD script
            Add all atoms to CPMD script
        """
        if self.qmatoms.empty:
            raise MiMiCPyError("No QM atoms have been selected")
        
        self.mpt.close() # clear all_data from memory, in case gc doesn't work soon enough
        ####
        #Write ndx, tpr file for gromacs
        #output only ndx if mdp is None
        ####
        
        inp = cpmd.Input()
        inp._ndx = self.qmatoms.index.to_list()
        
        # write index file
        if ndx != None:
            _global.host.write(_qmhelper.index(self.qmatoms.index, 'QMatoms'), ndx)
            _global.logger.write('info', f"Wrote Gromacs index file to {ndx}..")
            
        # sort by link column first, then element symbol
        # ensures that all link atoms are last, and all elements are bunched together
        # index is also reset, for getOverlaps_Atoms()
        sorted_qm = self.qmatoms.sort_values(by=['link', 'element']).reset_index()
        
        inp=cpmd.Input()
        
        inp.mimic = cpmd.Section()
        inp.mimic.paths = "1\n---" #path will be set by simulate.MiMiC.run()
        inp.mimic.box = '  '.join([str(s/bohr_rad) for s in self._mm_box])
        
        inp = _qmhelper.getOverlaps_Atoms(sorted_qm, inp) # get the overlaps and add atoms section of inp
        
        inp.mimic.long_range__coupling = ''
        inp.mimic.fragment__sorting__atom_wise__update = 100
        inp.mimic.cutoff__distance = 20.0
        inp.mimic.multipole__order = 3
        
        q = sum(self.qmatoms['charge'])
        if not round(q, 2).is_integer(): 
            _global.logger.write('warning', (f'Total charge of QM region (={q}) not an integer up to 2 decimal places.'
                                            ' \nRounding to integer anyways..'))
        
        inp.system.charge = round(q) # system section already created in getQverlap_Atoms()
        
        # set cpmd section
        if not inp.checkSection('cpmd'): inp.cpmd = cpmd.Section()
        inp.cpmd.mimic = ''
        inp.cpmd.parallel__constraints = ''
        
        inp.dft = cpmd.Section()
        inp.dft.functional__blyp = '' # TO DO: update to new XC_DRIVER code
        
        # default values
        inp.cpmd.maxsteps = 1000
        inp.cpmd.timestep = 0.001
        
        if inp_file is None:
            _global.logger.write('info', "Created CPMD input script..")
        else:
            _global.host.write(str(inp), inp_file)
            _global.logger.write('info', f"Wrote CPMD input script to {inp_file}..")
        
        return inp