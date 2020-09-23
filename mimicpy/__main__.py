#!/usr/bin/env python

import mimicpy
import sys
import time
import itertools
import threading
import os


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


def checkargs(args):
    ret = 0
    for k, v in args.items():
        if v is None:
            print(f"\n-{k} option should be specified!")
            ret += 1

    if ret:
        print(f'\nError in arguments passed. Exiting..\n')
        sys.exit(0)

def prepqm(qm, ndx, cpmd, cpmd_template, mdp):
    command = input('Please enter the selection below. Type help to get help message.\n>  ')
    command = command.split()
    prefix = command[0].lower()
    rest = " ".join(command[1:])

    if ndx == None: ndx = 'index.ndx'
    if cpmd == None: cpmd = 'cpmd.inp'

    if prefix == 'q'or prefix == 'quit':
        try:
            qm.getInp(cpmd_template, mdp, ndx, cpmd)
        except mimicpy.utils.errors.MiMiCPyError as e:
            print(e)
        return False
    elif prefix == 'add':
        try:
            qm.add(rest)
        except (mimicpy.utils.errors.SelectionError, mimicpy.utils.errors.MiMiCPyError) as e:
            print(e, end='\n\n')
        return True
    elif prefix == 'del':
        try:
            qm.delete(rest)
        except (mimicpy.utils.errors.SelectionError, mimicpy.utils.errors.MiMiCPyError) as e:
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
                    k, v = d.split()
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
        try:
            mpt_obj = mimicpy.parsers.MPT.fromTop(topol, nonstd_atm, mode='w')
        except (mimicpy.utils.errors.ParserError, FileNotFoundError) as e:
            print(e)
            return False
        mpt_obj.write(mpt)
    else:
        try:
            return mimicpy.parsers.MPT.fromTop(topol, nonstd_atm, mode='r')
        except (mimicpy.utils.errors.ParserError, FileNotFoundError) as e:
            print(e)
            return False


import argparse

def main():
    print('\n \t                ***** MiMiCPy CLI *****                  ')
    print('\n \t For more information type mimicpy [subcommand] --help \n')

    parser = argparse.ArgumentParser(prog='mimicpy')
    subparsers = parser.add_subparsers(title='valid subcommands',
                                       metavar='')  # Turns off list of subcommands

    parser_getmpt = subparsers.add_parser('getmpt',
                                          help='create MiMiCPy topology from GROMACS topology')
    getmpt_input_arguments = parser_getmpt.add_argument_group('options to specify input files')
    getmpt_input_arguments.add_argument('-top',  # TODO: Give better names for command line flags
                                        required=True,  # TODO: Make nicer indentations
                                        help='GROMACS topology file',
                                        metavar='topol.top')
    getmpt_output_arguments = parser_getmpt.add_argument_group('options to specify output files')
    getmpt_output_arguments.add_argument('-mpt',
                                         default='mimic.mpt',
                                         help='MiMiCPy topology file',
                                         metavar='mimic.mpt')
    parser_getmpt.set_defaults(func=getmpt)

    parser_prepqm = subparsers.add_parser('prepqm',
                                          help='create CPMD input and GROMACS index files')
    prepqm_input_arguments = parser_prepqm.add_argument_group('options to specify input files')
    prepqm_input_arguments.add_argument('-mpt',
                                        required=True,
                                        help='MiMiCPy mpt file',  # TODO: Print more precise help messages
                                        metavar='mimic.mpt')
    prepqm_input_arguments.add_argument('-gro', # TODO: Think of support for other formats
                                        required=True,
                                        help='GROMACS gro file',
                                        metavar='coords.gro')
    prepqm_input_arguments.add_argument('-inp',
                                        required=False,
                                        help='CPMD input script',
                                        metavar='cpmd.inp')
    prepqm_input_arguments.add_argument('-mdp',
                                        required=False,
                                        help='GROMACS input script',
                                        metavar='nve.mdp')
    prepqm_output_arguments = parser_prepqm.add_argument_group('options to specify output files')
    prepqm_output_arguments.add_argument('-out',
                                         default='mimic.inp',
                                         help='CPMD input script',
                                         metavar='mimic.inp')
    prepqm_output_arguments.add_argument('-ndx',
                                         default='mimic.ndx',
                                         help='GROMACS index file',
                                         metavar='index.ndx')
    parser_prepqm.set_defaults(func=prepqm)


    args = parser.parse_args()
    print(args)

if __name__ == '__main__':
    main()
