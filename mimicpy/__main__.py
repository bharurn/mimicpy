#!/usr/bin/env python

import argparse
import pandas as pd
import sys
import time
import itertools
import threading
import mimicpy


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


def prepqm(args):

    from mimicpy import Preparation, DefaultSelector

    def selection_help():
        print("\nvalid subcommands:\n\n"
              "    o add       <selection>    add selected atoms to QM region\n"
              "    o add-link  <selection>    add selected atoms to QM region as link atoms\n"
              "    o delete    <selection>    delete selected atoms from QM region\n"
              "    o clear                    clear all atoms from QM region\n"
              "    o view      <file-name>    print current QM region to console or file\n"
              "    o quit                     create CPMD input and GROMACS index files from selected atoms and quit\n"
              "    o help                     show this help message\n\n"
              "For more information on the selection langauge please refer to the docs.\n")

    def view(file_name=None):
        if prep.qm_atoms.empty:
            print("No QM atoms have been selected")
        else:
            # display whole dataframe
            with pd.option_context('display.max_rows', None, 'display.max_columns', None):
                    qmatoms_str = str(prep.qm_atoms)
                    
            if file_name:
                with open(file_name, 'w') as f:
                    f.write(qmatoms_str)
                print("Wrote list of QM atoms to {}".format(file_name))
            else:
                print(qmatoms_str)
                
    mpt = get_nsa_mpt(args)
    
    print('')
    loader = Loader('**Reading coordinates**')
    
    try:
        selector = DefaultSelector(mpt, args.coords)
        prep = Preparation(selector)
    except FileNotFoundError as e:
        print('\n\nError: Cannot find file {}! Exiting..\n'.format(e.filename))
        loader.close(halt=True)
        sys.exit()
    except mimicpy.utils.errors.ParserError as e:
        print(e)
        loader.close(halt=True)
        sys.exit()
        
    loader.close()
    
    dispatch = {'add':prep.add,
                'add-link': lambda selection: prep.add(selection, True),
                'delete':prep.delete,
                'clear': prep.clear,
                'view':view,
                'help':selection_help}
    print('\nPlease enter selection below. For more information type help')
    
    import readline

    readline.parse_and_bind('set editing-mode vi') # handle command history
        
    while True:
        user_input = input('> ')
        user_input = user_input.split()
        try:
            command = user_input[0].lower()
        except IndexError: # handle empty commands
            continue
        
        if command in ['quit', 'q']:
            try:
                if not prep.qm_atoms.empty:
                    prep.get_mimic_input(args.inp, args.mdp, args.ndx, args.out)
                break
            except mimicpy.utils.errors.SelectionError as error:
                print(error)
        selection = ' '.join(user_input[1:])
        try:
            dispatch[command](selection)
        except mimicpy.utils.errors.SelectionError as e:  # Selection errors
            print(e)
        except TypeError:  # include functions without argument
            dispatch[command]()
        except KeyError: # invalid commands
            print("Invalid command! Please try again. Type help for more information.")


def get_nsa_mpt(args):
    nsa_dct = {}
    if args.nsa:
        if args.top.split('.')[-1] == 'mpt':
            print("Non-standard atomtype file ignored as .mpt file was passed")
        else:
            print("\n**Reading non-standard atomtypes file**\n")
            from os.path import isfile
            if not isfile(args.nsa):
                print("Error: Cannot find file {}! Exiting..\n".format(args.nsa))
                sys.exit()
            else:
                with open(args.nsa, 'r') as f:
                    for i, line in enumerate(f.readlines()):
                        splt = line.split()
                        if len(splt) < 2:
                            print("Line {} in nonstandard atomtypes file {} not in 2-column format!\n".format(i+1, args.nsa))
                            sys.exit()
                        elif len(splt) > 2:
                            print("Line {} in nonstandard atomtypes file {} has more than 2-columns. Using first two values only.\n".format(i+1, args.nsa))
                        
                        nsa_dct[splt[0]] = splt[1]
            
    if nsa_dct != {}:
        print("The following atomypes were read from {}:\n".format(args.nsa))
        print("+---------------------+")
        print("|{:^11}|{:^9}|".format("Atom Type", "Element"))
        print("+---------------------+")
        for k,v in nsa_dct.items():
            print("|{:^10} | {:^8}|".format(k, v))
            print("+---------------------+")
    
    print("\n**Readining topology**\n")
    
    try:
        return mimicpy.Mpt.from_file(args.top, mode='w', nonstandard_atomtypes=nsa_dct)
    except FileNotFoundError as e:
        print('\n\nError: Cannot find file {}! Exiting..\n'.format(e.filename))
        sys.exit()

def getmpt(args):
    get_nsa_mpt(args).write(args.mpt)


def main():
    print('\n \t                ***** MiMiCPy *****                  ')
    print('\n \t For more information type mimicpy [subcommand] --help \n')

    parser = argparse.ArgumentParser(prog='mimicpy')
    subparsers = parser.add_subparsers(title='valid subcommands',
                                       metavar='')  # Turns off list of subcommands

    parser_getmpt = subparsers.add_parser('getmpt',
                                          help='create MiMiCPy topology from GROMACS topology')
    getmpt_input = parser_getmpt.add_argument_group('options to specify input files')
    getmpt_input.add_argument('-top',
                              required=True,
                              help='GROMACS topology file',
                              metavar='topol.top')
    getmpt_input.add_argument('-nsa',
                              required=False,
                              help='list of non-standard atomtypes',
                              metavar='')
    getmpt_output = parser_getmpt.add_argument_group('options to specify output files')
    getmpt_output.add_argument('-mpt',
                               default='topol.mpt',
                               help='MiMiCPy topology file',
                               metavar='topol.mpt')
    parser_getmpt.set_defaults(func=getmpt)

    parser_prepqm = subparsers.add_parser('prepqm',
                                          help='create CPMD MiMiC input and Gromacs index files')
    prepqm_input = parser_prepqm.add_argument_group('options to specify input files')
    prepqm_input.add_argument('-top',
                              required=True,
                              help='Topology file',
                              metavar='[.top/.mpt]')
    prepqm_input.add_argument('-coords',
                              required=True,
                              help='Coordinate file',
                              metavar='[.gro/.pdb]')
    prepqm_output = parser_prepqm.add_argument_group('options to specify output files')
    prepqm_output.add_argument('-out',
                               default='cpmd.inp',
                               help='CPMD script for MiMiC run',
                               metavar='[.inp] (cpmd.inp)')
    prepqm_output.add_argument('-ndx',
                               default='index.ndx',
                               help='Gromacs index file',
                               metavar='[.ndx] (index.ndx)')
    prepqm_others = parser_prepqm.add_argument_group('other options')
    prepqm_others.add_argument('-nsa',
                              required=False,
                              help='list of non-standard atomtypes in 2-column format',
                              metavar='[.txt/.dat]')
    prepqm_others.add_argument('-inp',
                              required=False,
                              help='CPMD template input script',
                              metavar='[.inp]')
    prepqm_others.add_argument('-mdp',
                              required=False,
                              help='Gromacs template MDP script',
                              metavar='[.mdp]')
    parser_prepqm.set_defaults(func=prepqm)

    args = parser.parse_args()
    if vars(args) == {}:
        sys.exit()
    subcommand = args.func.__name__
    print('=====> Running {} <=====\n'.format(subcommand))
    args.func(args)
    print('\n=====> Done <=====\n'.format(subcommand))
    
if __name__ == '__main__':
    main()
