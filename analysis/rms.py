from .plotbox import plotBoxGen
import pandas as pd

@plotBoxGen('RMSD', 'Step', 'Time')
def d(top, trr, selections):
    try:
        import MDAnalysis as mda
        from MDAnalysis.analysis import align
    except:
        raise ImportError("Please install MDAnalysis python package to parse trajectory files.")
        
    u = mda.Universe(top)

    ref = mda.Universe(top)
    
    rmsd = []
    
    for l in trr:
        u.trajectory =  mda.Universe(l).trajectory
    
        rmsd.extend([(i.data['step'], i.data['step']*i.dt/10000, *[align.alignto(u.select_atoms(s), ref.select_atoms(s))[1]\
                      for s in selections]) for i in u.trajectory])
        u.trajectory.close()
    
    cols = [' '.join([w.title() if w.islower() else w for w in s.split()]) for s in selections]
    df = pd.DataFrame(rmsd, columns=('Step', 'Time', *cols))
    return df.set_index(['Step', 'Time'])

@plotBoxGen('RMSF', 'Step', 'Time')
def f(top, trr, selection, selections):
    
    try:
        import MDAnalysis as mda
        from MDAnalysis.analysis import align
    except:
        raise ImportError("Please install MDAnalysis python package to parse trajectory files.")
    
    u = mda.Universe(top, trr)
    align.AlignTraj(u, u, select=selection, in_memory=True).run()
    