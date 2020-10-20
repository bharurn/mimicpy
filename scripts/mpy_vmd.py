import sys
import mimicpy

#TODO: error handling

if __name__ == '__main__':
    mpt = sys.argv[1]
    cpmd = sys.argv[2]
    ndx = sys.argv[3]
    
    qm = mimicpy.Prepare(mimicpy.selector.VMD(mpt, tcl_vmd_params=sys.argv[4:]))
    qm.add() # passed selection doesn't matter
    qm.get_mimic_input(ndx_file=ndx, cpmd_file=cpmd)