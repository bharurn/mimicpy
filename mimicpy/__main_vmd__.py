#!/usr/bin/env python

import sys
import mimicpy

def main():
    if len(sys.argv) < 22:
        print("Not enough arguments passed. Exiting..\n")
        sys.exit(1)
    
    top = sys.argv[1]
    inp = None if sys.argv[2] == 'None' else sys.argv[2]
    mdp = None if sys.argv[3] == 'None' else sys.argv[3]
    ndx = sys.argv[4]
    out = sys.argv[5]
    #sys.argv[6:] should have all selection info from VMD
    
    try:
        qm = mimicpy.Prepare(mimicpy.selector.VMD(top, tcl_vmd_params=sys.argv[6:]))
    except FileNotFoundError as e:
        print('\n\nError: Cannot find file {}! Exiting..\n'.format(e.filename))
        sys.exit(1)
    except (mimicpy.utils.errors.ParserError, mimicpy.utils.errors.MiMiCPyError) as e:
        print(e)
        sys.exit(1)
    
    try:
        qm.add() # passed selection doesn't matter
    except mimicpy.utils.errors.MiMiCPyError as e:
        print(e)
        sys.exit(1)
    
    try:
        qm.get_mimic_input(inp, mdp, ndx, out)
    except FileNotFoundError as e:
        print('\n\nError: Cannot find file {}! Exiting..\n'.format(e.filename))
        sys.exit(1)
    except mimicpy.utils.errors.SelectionError as e:
        print(e)

if __name__ == '__main__':
    main()