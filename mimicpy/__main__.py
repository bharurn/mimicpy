#!/usr/bin/env python

import logging
import mimicpy
import sys
import time
import itertools
import threading


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


#def prepqm(qm, ndx, cpmd, cpmd_template, mdp):
def prepqm(args):
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
              "For more information on the selection langauge please refer to the docs.\n") # What docs :P
    else:
        print("Command not understood! Type help to get a list of accepted commands.\n")
        return True

def prepqm(args):
    from mimicpy.core.prepare import Preparation
    help = lambda _ : print("\nvalid subcommands:\n\n"
                            "    add <selection>    add selected atoms to QM region\n"
                            "    delete <selection>    delete selected atoms from QM region\n"
                            "    clear              clear all atoms from QM region\n"
                            "    view               print current QM region to console\n"
                            "    quit               create CPMD input and GROMACS index files from selected atoms and quit\n"
                            "    help               show this help message\n\n"
                            "For more information on the selection langauge please refer to the docs.\n") # What docs :P
    prep = Preparation(args.mpt, args.gro)
    dispatch = {'add':prep.add,
                'delete':prep.delete,
                'clear':prep.clear,
                'view':prep.view}
    print('Please enter selection below. For more information type help')
    while True:
        user_input = input('> ')
        user_input = user_input.split()
        command = user_input[0].lower()
        if command == 'quit':
            break
        selection = ' '.join(user_input[1:])
        try:
            dispatch[command](selection)
        except:  # TODO: Write appropriate error message
            pass

def getmpt(args):
    # TODO: Write reader for non-standard atomtypes (io)
    # TODO: Interface non-standard atomtype file (in top.__read)
    from mimicpy.io.mpt import Mpt
    Mpt.from_file(args.top, mode='w', nonstandard_atomtypes=None).write(args.mpt)


import argparse

def main():
    logging.info('Invoked new mimicpy command')

    print('\n \t                ***** MiMiCPy CLI *****                  ')
    print('\n \t For more information type mimicpy [subcommand] --help \n')

    parser = argparse.ArgumentParser(prog='mimicpy')
    subparsers = parser.add_subparsers(title='valid subcommands',
                                       metavar='')  # Turns off list of subcommands

    parser_getmpt = subparsers.add_parser('getmpt',
                                          help='create MiMiCPy topology from GROMACS topology')
    getmpt_input = parser_getmpt.add_argument_group('options to specify input files')
    getmpt_input.add_argument('-top',
                              required=True,  # TODO: Make nicer indentations
                              help='GROMACS topology file',
                              metavar='topol.top')
    getmpt_input.add_argument('-nsa',  # TODO: Give better names for command line flags
                              required=False,
                              help='non-standard atomtypes',
                              metavar='')
    getmpt_output = parser_getmpt.add_argument_group('options to specify output files')
    getmpt_output.add_argument('-mpt',
                               default='mimic.mpt',
                               help='MiMiCPy topology file',
                               metavar='mimic.mpt')
    parser_getmpt.set_defaults(func=getmpt)

    parser_prepqm = subparsers.add_parser('prepqm',
                                          help='create CPMD input and GROMACS index files')
    prepqm_input = parser_prepqm.add_argument_group('options to specify input files')
    prepqm_input.add_argument('-mpt',
                              required=True,
                              help='MiMiCPy mpt file',  # TODO: Print more precise help messages
                              metavar='mimic.mpt')
    prepqm_input.add_argument('-gro', # TODO: Think of support for other formats
                              required=True,
                              help='GROMACS gro file',
                              metavar='coords.gro')
    prepqm_input.add_argument('-inp',
                              required=False,
                              help='CPMD input script',
                              metavar='cpmd.inp')
    prepqm_input.add_argument('-mdp',
                              required=False,
                              help='GROMACS input script',
                              metavar='nve.mdp')
    prepqm_output = parser_prepqm.add_argument_group('options to specify output files')
    prepqm_output.add_argument('-out',
                               default='mimic.inp',
                               help='CPMD input script',
                               metavar='mimic.inp')
    prepqm_output.add_argument('-ndx',
                               default='mimic.ndx',
                               help='GROMACS index file',
                               metavar='index.ndx')
    parser_prepqm.set_defaults(func=prepqm)

    args = parser.parse_args()
    subcommand = args.func.__name__
    logging.info(f'Started mimicpy {subcommand}.')
    args.func(args)  # TODO: Put annimation back
    logging.info(f'Finished mimicpy {subcommand}.')
    print(f'MiMiCPy {subcommand} finished successfully')

if __name__ == '__main__':
    main()
