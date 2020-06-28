import sys
import mimicpy

#TO DO: error handling

mimicpy.setLogger(1)

if __name__ == '__main__':
    mpt = sys.argv[1]
    gro = sys.argv[2]
    cpmd = sys.argv[3]
    ndx = sys.argv[4]
    
    qm = mimicpy.Prepare(mpt, None, selector=mimicpy.selector.VMD(tcl_vmd_params=sys.argv[5:]))
    qm.add() # passed selection doesn't matter
    qm.getInp(ndx, cpmd)