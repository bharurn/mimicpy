#import mimicpy
import argparse
import sys

class MiMiCPyParser(argparse.ArgumentParser):
    def __init__(self):
        super().__init__(usage='MiMiCPy CLI Program')
        
    def error(self, message):
        print("NOOOO!!!")

def interact():
    command = input('Please enter the selection below. Type help to get help message.\n>  ')
    command = command.split()
    prefix = command[0].lower()
    if prefix == 'q'or prefix == 'quit':
        return False
    elif prefix == 'add':
        return True
    elif prefix == 'del':
        return True
    elif prefix == 'clear':
        return True
    elif prefix == 'help' or prefix == 'h':
        print("Following commands are understood:\n\n"
              "- add <selection>\n\tAdd atoms (given by the selection) to the QM region\n\n"
              "- del <selection>\n\tDelete atoms (given by selection) from the QM region\n\n"
              "- clear\n\tClear all atoms from the QM region to start over\n\n"
              "- quit or q\n\tGenerate the CPMD input and Gromacs index files from selected atoms, and quit.\n\n"
              "- help or h\n\tShow this help message.\n\n"
              "A keyboard interrupt will cause an immediate termination. "
              "Also, for more information on the selection langauge please refer to the docs.\n")
    else:
        print("Command not understood! Type help to get a list of accepted commands.\n")
        return True

if __name__ == '__main__':
    
    parser = MiMiCPyParser()
    parser.add_argument("subprogram", help="Subprogram to be called", type=str)
    
    args = parser.parse_args()
    
    if args.subprogram != 'prepqm':
        sys.exit(0)
    while True:
        try:
            ret = interact()
        except (KeyboardInterrupt, EOFError):
            print("\nExiting without saving..")
            break
        if ret == False: break