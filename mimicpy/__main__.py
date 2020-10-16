#!/usr/bin/env python

import logging
import argparse
import sys
import time
import itertools
import threading
import mimicpy  # TODO: Adjust __init__ files (leave now for logging config)


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


def prepqm(args):  # TODO: accept non-standard atomtypes - ignore if mpt is passsed

    from mimicpy.core.prepare import Preparation
    from mimicpy.core.selector import GroSelector

    def selection_help():  # TODO: Find a better way or place for help
        print("\nvalid subcommands:\n\n"
              "    o add <selection>    add selected atoms to QM region\n"
              "    o delete <selection> delete selected atoms from QM region\n"
              "    o clear              clear all atoms from QM region\n"
              "    o view               print current QM region to console\n"
              "    o quit               create CPMD input and GROMACS index files from selected atoms and quit\n"
              "    o help               show this help message\n\n"
              "For more information on the selection langauge please refer to the docs.\n") # What docs :P

    def view():
        print(prep.get_qm_atoms())  # TODO: Format printing

    selector = GroSelector(args.mpt, args.gro)
    prep = Preparation(selector)
    dispatch = {'add':prep.add,  # TODO: Add as link atom
                'delete':prep.delete,
                'clear': prep.clear,
                'view':view,
                'help':selection_help}
    print('Please enter selection below. For more information type help')
    while True:
        user_input = input('> ')
        user_input = user_input.split()
        command = user_input[0].lower()
        if command in ['quit', 'q']:
            try:
                prep.get_mimic_input(args.inp, args.mdp, args.ndx, args.out)
                break
            except mimicpy.utils.errors.SelectionError as error:
                print(error)
        selection = ' '.join(user_input[1:])
        try:
            dispatch[command](selection)  # TODO: Give feedback
        except TypeError:  # TODO: Write appropriate error message
            dispatch[command]()
        except:
            pass


def getmpt(args):
    # TODO: Write reader for non-standard atomtypes (io)
    # TODO: Interface non-standard atomtype file (in top.__read)
    from mimicpy.io.mpt import Mpt
    Mpt.from_file(args.top, mode='w', nonstandard_atomtypes=None).write(args.mpt)


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
    if vars(args) == {}:
        sys.exit()
    subcommand = args.func.__name__
    logging.info('Started mimicpy %s.', subcommand)
    args.func(args)  # TODO: Put annimation back
    logging.info('Finished mimicpy %s.', subcommand)
    print(f'MiMiCPy {subcommand} finished successfully')

if __name__ == '__main__':
    main()
