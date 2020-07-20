#!/usr/bin/python

import mimicpy
import sys
import time
import itertools
import threading
import os

mimicpy.setLogger(1)

class Loader:
    def __init__(self, message):
        self.done = False
        self.message = message
        t = threading.Thread(target=self.__animate)
        t.start()

    def __animate(self):
        for c in itertools.cycle(['|', '/', '-', '\\']):
            if self.done:
                break
            sys.stdout.write(f'\r{self.message}  ' + c)
            sys.stdout.flush()
            time.sleep(0.1)
    
    def close(self, halt=False):
        self.done = True
        if not halt: print('Done')
        
class MiMiCPyParser():
    
    subprograms = {'getmpt': 'Generate MiMiCPy topology from Gromacs topology', 
               'prepqm': 'Prepare QM region by generating the CPMD input and Gromacs index files',
               'cpmd2coords': 'Convert the MIMIC/ATOMS section of CPMD input file to a pdb/gro file'}
    
    subprogram_help = {'getmpt': "\n-top\n\tInput Gromacs topology file\n"
                       "-mpt\n\tOutput MiMiCPy topology file\n"
                       "\nOptional arguments:\n"
                       "\n-nonstd\n\tNon-standard residues list as text file\n"
                       "-pdb\n\tpdb file to read non-standard residue elements from"
                       "\n\tIt must have the element column filled for those residues!"
                       "\n-itp\n\titp file to read non-standard residue atom types from"
                       "\n-reslist\n\tList of nonstandard residues to read from pdb/itp files\n"
                       "\nThe last three options must be given together\n",
                'prepqm': "\nInput options:\n\n"
                      "-top\n\tGromacs topology file\n"
                      "-mpt\n\tMiMiCPy topology file\n"
                      "\nOnly one of the above need to be specified\n\n"
                      "-gro\n\tGromacs coordinate file\n"
                      "\nOutput options:\n\n"
                      "-cpmd\n\tCPMD input script name, deafult: cpmd.inp\n"
                      "-ndx\n\tGromacs index file name, default: index.ndx\n"
                      "\nOptional arguments for -top option:\n"
                      "\n-nonstd\n\tNon-standard residues list as text file\n"
                       "-pdb\n\tpdb file to read non-standard residue elements from"
                       "\n\tIt must have the element column filled for those residues!"
                       "\n-itp\n\titp file to read non-standard residue atom types from"
                       "\n-reslist\n\tList of nonstandard residues to read from pdb/itp files\n"
                       "\nThe last three options must be given together\n",
            'cpmd2coords': "\n-cpmd\n\tCPMD input script with MIMIC and ATOMS sections\n"
                            "\n-mpt\n\tMiMiCPy topology\n"
                            "\n-coords\n\tName of output coordinate file, defaults to conf.pdb\n"}
    
    def __init__(self):
        self.program = sys.argv[0]
        
        if len(sys.argv) == 1:
            self.help()
        else:
            self.subprogram = sys.argv[1]
        
        if self.subprogram not in self.subprograms.keys():
            self.help('Invalid subprogram!')
        
        args = [i for i in sys.argv[2::2]]
        vals = [i for i in sys.argv[3::2]]
        
        if args == [] or args == ['-help'] or args == ['-h']:
            print(f"mimicpy {self.subprogram}: {self.subprograms[self.subprogram]}")
            print(self.subprogram_help[self.subprogram])
            sys.exit(0)
        else:
            print(f"MiMiCPy CLI subprogram requested: {self.subprogram}\n")
        
        for k,v in zip(args, vals):
            setattr(self, k[1:], v)
    
    def __getattr__(self, key):
        if key not in self.__dict__:
            return None
        
    def help(self, message=f'\tUsing MiMiCPy version {mimicpy.__version__}'):
        
        print(f"\n{message}\n\n"
              "List of available subprograms:\n")
        
        for k,v in self.subprograms.items():
            print(f"\no {k}:\n\t{v}")
              
        print(f"\nFormat:\n\tmimicpy <subprogram> <options for subprogram>\n")
        sys.exit(0)

def checkargs(args):
    ret = 0
    for k,v in args.items():
        if v is None:
            print(f"\n-{k} option should be specified!")
            ret += 1
    
    if ret:
        print(f'\nError in arguments passed. Exiting..\n')
        sys.exit(0)

def prepqm(qm, ndx, cpmd):
    command = input('Please enter the selection below. Type help to get help message.\n>  ')
    command = command.split()
    prefix = command[0].lower()
    rest = " ".join(command[1:])
    if ndx == None: ndx = 'index.ndx'
    if cpmd == None: cpmd = 'cpmd.inp'
    
    if prefix == 'q'or prefix == 'quit':
        try:
            qm.getInp(ndx, cpmd)
        except mimicpy.utils.errors.MiMiCPyError:
            pass
        return False
    elif prefix == 'add':
        try:
            qm.add(rest)
        except mimicpy.utils.errors.SelectionError as e:
            print(e, end='\n\n')
        return True
    elif prefix == 'del':
        try:
            qm.delete(rest)
        except mimicpy.utils.errors.SelectionError as e:
            print(e, end='\n\n')
        return True
    elif prefix == 'clear':
        qm.clear()
        return True
    elif prefix == 'help' or prefix == 'h':
        print("Following commands are understood:\n\n"
              "o add <selection>\n\tAdd atoms (given by the selection) to the QM region\n\n"
              "o del <selection>\n\tDelete atoms (given by selection) from the QM region\n\n"
              "o clear\n\tClear all atoms from the QM region to start over\n\n"
              "o quit or q\n\tGenerate the CPMD input and Gromacs index files from selected atoms, and quit.\n\n"
              "o help or h\n\tShow this help message.\n\n"
              "A keyboard interrupt will cause an immediate termination. "
              "For more information on the selection langauge please refer to the docs.\n")
    else:
        print("Command not understood! Type help to get a list of accepted commands.\n")
        return True

def getmpt(topol, mpt, nonstd, reslist, pdb, itp):
    print("Parsing topology..")
    nonstd_atm = {}
    
    if nonstd:
        if not os.path.isfile(nonstd): raise FileNotFoundError(f"{nonstd} not found!")
        with open(nonstd, 'r') as f:
            for d in f.read().splitlines():
                try:
                    k,v = d.split()
                except ValueError:
                    print("Nonstandard residue file not in the correct format!\n")
                    sys.exit(0)
                nonstd_atm[k] = v
        
    elif reslist or pdb or itp:
        checkargs({'reslist': reslist, 'pdb': pdb, 'itp': itp})
        if not os.path.isfile(reslist): raise FileNotFoundError(f"{reslist} not found!")
        resname = []
        with open(reslist, 'r') as f:
            resname = f.read().splitlines()
            loader = Loader('Reading coordinates')
            nonstd_atm = mimicpy.parsers.top.nonStdTypes(pdb, itp, buff=1000, *resname)
            loader.close()
            
    if mpt:
        mpt_obj = mimicpy.parsers.MPT.fromTop(topol, nonstd_atm, mode='w')
        mpt_obj.write(mpt)
    else:
        return mimicpy.parsers.MPT.fromTop(topol, nonstd_atm, mode='r')

def main():
    print('\n\t   ****MiMiCPy CLI****\n')
    
    parser = MiMiCPyParser()
    
    if parser.subprogram == 'getmpt':
        checkargs({'mpt': parser.mpt, 'top': parser.top})
        try:
            getmpt(parser.top, parser.mpt, parser.nonstd, parser.reslist, parser.pdb, parser.itp)
        except Exception as e:
            print(e)
            print("Exiting with error..\n")
            sys.exit(0)
        print("MPT succesfully prepared..\n")
        
    elif parser.subprogram == 'prepqm':
        
        if parser.mpt is not None:
            mpt = parser.mpt
        elif parser.top is not None:
            mpt = getmpt(parser.top, None, parser.nonstd, parser.reslist, parser.pdb, parser.itp)
        else:
            print("Either -top or -mpt must be specified as an input file\n")
            sys.exit(0)
            
        checkargs({'gro': parser.gro})
            
        print("Preparing QM region..")
        
        loader = Loader('Reading topology and coordinates')
        
        try:
            qm = mimicpy.Prepare(mpt, parser.gro)
        except Exception as e:
            print(e)
            print("Exiting with error..\n")
            loader.close(halt=True)
            sys.exit(0)
        
        loader.close()
        
        while True:
            try:
                ret = prepqm(qm, parser.ndx, parser.cpmd)
            except (KeyboardInterrupt, EOFError):
                print("\nExiting without saving..")
                ret = False
                
            if ret == False:
                print('')
                break
            
    elif parser.subprogram == 'cpmd2coords':
        
        checkargs({'mpt': parser.mpt, 'cpmd': parser.cpmd})
        cpmd = mimicpy.scripts.cpmd.Input.fromFile(parser.cpmd)
        
        loader = Loader('Reading topology')
        mpt = mimicpy.parsers.MPT.fromFile(parser.mpt)
        loader.close()
        
        if parser.coords == None: parser.coords = 'conf.pdb'
        cpmd.toCoords(mpt, parser.coords)
        
        print(f"Wrote output to {parser.coords}")
