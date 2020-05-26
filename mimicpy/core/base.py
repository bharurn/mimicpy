# This is part of MiMiCPy

"""

This module contains the BaseHandle class which is inherited by all handles

"""

import re
import yaml
from .._global import _Global as _global
from ..utils.errors import GromacsError, MiMiCPyError, EnvNotSetError, defaultHook
from ..utils.logger import Logger, LogString
import sys

class BaseHandle:
    """
    Contains functions/variables common to all handles:
        prepare.MM, prepare.QM, simulate.MD, simulate.MiMiC
    All handles inherit from this class
    
    """
    
    def __init__(self, status=None):
        """Init status and log string"""
        if not status:
            status = {'prepMM': '', 'prepQM': '', 'run': ['']}
        self._status = status # init _status
        self.log = LogString() # log string of standard ouput from all gmx commands
        # init logger with gmx log string, and notes redirected to stdout
        self.logger = Logger(log=self.log, notes=sys.stdout)
        self.current_cmd = 'gmx'
    
    def getcurrent(self, ext, level=False, exp=True):
        """Find the latest file with extnestion ext using _status"""
        
        _dir = _global.host.pwd()+'/' # get cwd
        
        if ext == 'top' or ext == 'mpt':
            # if top/mpt go directly to prepMM folder
            return _dir+BaseHandle._getFile(self._status['prepMM'], ext)
        elif ext == 'mimic-tpr':
            # if mimic tpr is asked for go directly to prepQM folder
            return _dir+BaseHandle._getFile(self._status['prepQM'], ext)
        
        # if topology is not asked for, then it is run files (trr,cpt,etc.)
        # get the list of folders from run key, also tack on the prepMM and prepQM folder
        # just in case pdb, etc. files were asked for
        run =  [self._status['prepMM'], self._status['prepQM']] + self._status['run']
        
        for i,d in enumerate(run[::-1]): # loop in reverse, latests to earliest
            ret = BaseHandle._getFile(d, ext)
            
            if ret != None: # if we got a file, return it, else continue with next directory
                if level: return(i, _dir+ret)
                else: return _dir+ret
        
        # if nothing was found in any folder, raise exception
        if exp: raise FileNotFoundError(f"Cannot find file with extension {ext}")
    
    def __del__(self):
        """Deconstructor"""
        self.logger.close()
        self.toYaml()
        
    def getcurrentNone(self, file, ext, level=False, exp=True):
        """
        Convenience function, check if argument passed is None
        if it is searches for file in folders
        else it just return back that file
        """
        if file:
            return file
        else:
            return self.getcurrent(ext, level, exp)
    
    def gethistory(self, ext):
        _dir = _global.host.pwd()+'/' # get cwd
        
        lst = []
        run =  [self._status['prepMM'], self._status['prepQM']] + self._status['run']
        
        for i,d in enumerate(run[::-1]): # loop in reverse, latests to earliest
            ret = BaseHandle._getFile(d, ext)
            
            if ret != None: # if we got a file, return it, else continue with next directory
                lst.append(_dir+ret)
        
        return lst
    
    @staticmethod
    def _getFile(dirc, ext):
        """Finds file in dirc folder with extension ext"""
        
        # get list of all files in dirc with extension ext, ignore folders
        lst = _global.host.ls(dirc=dirc, file_eval=lambda a: True if a.endswith(ext) else False, dir_eval=lambda a: False)
        
        if ext == 'cpt': # gromacs saves multiple cpt files, take only the first one
            lst = [l for l in lst if '_prev' or '_step' not in l]
        
        if lst == []: return None # if nothing was found
        # more than one was found, raise exception
        elif len(lst) > 1: raise MiMiCPyError(f"More than one current file found with extension {ext}!"
                                        f"\nFiles found: {','.join(lst)}")
        if dirc.strip() != '': dirc += '/'
            
        return dirc+lst[0] # return with the directory
    
    def getStatus(self): return self._status
    
    @classmethod
    def continueFrom(cls, session, *args, **kwargs):
        """Transfer _status from session to new handle"""
        if not hasattr(session, '_status'):
            raise MiMiCPyError("No simulation status variables were found for session!")
        
        return cls(status=session._status, *args, **kwargs)
    
    def toYaml(self):
        """Save _status to yaml"""
        _global.logger.write('debug', f"Saving status to _status.yaml..")
        y = yaml.dump(self._status)
        _global.host.write(y, '_status.yaml')
    
    @classmethod
    def fromYaml(cls, *args, **kwargs):
       """Transfer _status from _status.yaml to new handle"""
       _global.logger.write('debug', f"Loading status from _status.yaml..")
       txt = _global.host.read('_status.yaml')
       _status = yaml.safe_load(txt)
       return cls(status=_status, *args, **kwargs)
   
    def setcurrent(self, dirc=None, key='run'):
        """Add new directory to _status dict"""
        
        if dirc == None: # if the argument is None add self.dir
            # meant for prepare classes
            dirc = self.dir
        _global.host.mkdir(dirc)
        if key == 'prepMM' or key == 'prepQM': # add the correct key
            self._status[key] = dirc
        elif dirc not in self._status['run']: # add to run list
            self._status['run'].append(dirc)
    
    def _gmxhook(self, cmd, text):
        """
        Function called by _global.host eveytime a gmx command is executed
        Used to parse gmx ouput and check for errors/notes/warnings
        """
        defaultHook(cmd, text) # call defaultHook to look for stupid errors
        
        # write log
        self.logger.write('log', f"==>Command Run: {cmd}\n")
        self.logger.write('log', text)
        
        self._gmxerrhdnl(cmd, text) # check for errors
        
        notes = BaseHandle._notes(text) # get notes/warnings
        
        if not notes.isspace() and notes.strip() != '': # write notes
            self.logger.write('notes', f"From {self.current_cmd}")
            self.logger.write('notes', notes)
    
    @staticmethod
    def _gmxerrhdnl(gmx_cmd, text, dont_raise=False):
        """Parse gmx output for errors and raise exception"""
        
        # checks for these words in output; should be expanded
        if 'Halting program' in text or 'Fatal error' in text or 'Error in user input' in text:
            
            pattern = re.compile("^-+$")
            
            msg = []
            start = False
            
            for line in text.splitlines()[::-1]: # read backwards
                if line.startswith('For more information and tips for troubleshooting') and not start:
                    start=True
                elif start:
                    if pattern.match(line):
                        break
                    else:
                        msg.append(line)
                    
                err = '\n'.join(msg[::-1])
                  
            if dont_raise:
                err += f"Error running Gromacs {gmx_cmd}!\n{err}"
            else:
                raise GromacsError(gmx_cmd, err) # raise exception
        
        if dont_raise:
            return err
        else:
            return None
        
    @staticmethod
    def _notes(text):
        """Parse gmx output for notes and warnings"""
        
        pattern = re.compile("^-+$")
        
        notes = ''
        start = False
        for i, line in enumerate(text.splitlines()):
            # parse from here
            if line.startswith('NOTE') or line.startswith('WARNING') or 'One or more water molecules can not be settled.' in line:
                notes += line + '\n'
                start = True
            elif start:
                if line.strip() == '' or line.startswith('WARNING') or line.startswith('NOTE') or pattern.match(line):
                    # end parsing
                    start = False
                else:
                    notes += line + '\n'
        
        return ''.join(filter(lambda x: not re.match(r'^\s*$', x), notes)) # remove blank lines
    
     
    def gmx(self, cmd, **kwargs):
        """
        Function to execute gmx commands in "pythonic" way
        e.g., to execute gmx pdb2gmx -f in.pdb -o out.gro -water spce -his
        call gmx('pdb2gmx', f='in.pdb', o='out.gro', water='spce', his='')
        """
        
        if _global.gmx is None or _global.gmx.strip() == '': # make sure gmx path is set
            raise EnvNotSetError('Gromacs executable', 'gmx')
            
        gmx_ex = _global.gmx.strip()
        
        if 'onlycmd' in kwargs: # if onlycmd is passed set it to true
            onlycmd = True
            del kwargs['onlycmd']
        else:
            onlycmd = False
        
        if 'dirc' in kwargs: # if dirc is passed set it
            dirc = kwargs['dirc']
            del kwargs['dirc']
        else:
            dirc = ''
        
        # if its mdrun, we need to run it in the background
        if cmd == 'mdrun': mdrun = True
        else: mdrun = False
        
        cmd = f"{gmx_ex} -quiet {cmd}"
        
        stdin = None
        
        # parse arguments to convert into gromacs command
        for k, v in kwargs.items():
            if k != 'stdin':
                cmd += f" -{k} {v}"
            elif k == 'stdin':
                stdin = v
                
        if onlycmd: # if onlycmd was set, don't run just return the cmd string
            # used mainly with Slurm
            return cmd
        
        _global.logger.write('debug', f"Running {cmd}..")
        
        self.current_gmx = cmd
        if mdrun: _global.host.runbg(cmd, hook=self._gmxhook, dirc=dirc) # runbg for mdrun
        else: _global.host.run(cmd, stdin=stdin, hook=self._gmxhook, dirc=dirc) # run for everything else
    
    def grompp(self, mdp, new, **kwargs):
        """
        Writes mdp to MDP file
        Finds latest gro/trr, and top
        Execute gmx grompp
        """
        new = new.lower().replace(' ', '_')
        # write mdp file into folder dirc
        if 'dirc' in kwargs:
            mdp_file = f'{kwargs["dirc"]}/{new}.mdp'
            # don't delete kwargs['dirc'], it should be passed to gmx
        else:
            mdp_file = f'{new}.mdp'
        _global.host.write(str(mdp), mdp_file)
        
        # a custom gro/trr file can be passed to grompp
        if 'gro' in kwargs:
            gro_file = kwargs['gro']
            del kwargs['gro']
        else:
            gro_file = None
            
        if 'trr' in kwargs:
            trr_file = kwargs['trr']
            del kwargs['trr']
        else:
            trr_file = None
        
        # if no custom file is passed, search for it in folders
        gro = self.getcurrentNone(gro_file, 'gro', level=True, exp=False)
        trr = self.getcurrentNone(trr_file, 'trr', level=True, exp=False)
        
        if gro == None and trr ==  None:
            raise MiMiCPyError("No coordinate data (gro/trr file) was found..")
        
        if trr != None: # if both gro and trr were found
            if gro[0] < trr[0]: # find the latest one, and run grompp accordingly
                self.gmx('grompp', f = f'{new}.mdp', c = gro[1],\
                     p = self.getcurrent('top'), o = f"{new}.tpr", **kwargs)
            else:
                self.gmx('grompp', f = f'{new}.mdp', c = gro[1],\
                     p = self.getcurrent('top'), t=trr[1], o = f"{new}.tpr", **kwargs)
        else: # otherwise, just use the gro file
            self.gmx('grompp', f = f'{new}.mdp', c = gro[1],\
                     p = self.getcurrent('top'), o = f"{new}.tpr", **kwargs)
            
    def cpmd(self, inp, out, onlycmd=False, dirc=''):
        """
        Function to execute gmx commands in "pythonic" way
        e.g., to execute cpmd cpmd.in path/to/pp > cpmd.oput
        call cpmd('cpmd.in', 'cpmd.out')
        Make sure cpmd_pp is set in _Global before that!
        """
        
        # check for env variables
        if _global.cpmd is None or _global.cpmd.strip() == '':
            raise EnvNotSetError('CPMD executable', 'cpmd')
        
        if _global.cpmd_pp is None or _global.cpmd_pp.strip() == '':
            raise EnvNotSetError('CPMD pseudopotential', 'cpmd_pp')
        
        cmd = f"{_global.cpmd} {inp} {_global.cpmd_pp} > {out}"
        
        if onlycmd: return cmd # return only cmd
        
        _global.logger.write('debug', cmd)
        
        _global.logger.write('debug', "Running {cmd}..")
        _global.host.runbg(cmd, dirc=dirc)