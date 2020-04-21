from pygmx.shell import remote, local

cmd = local.Local()

def set(work_dir, modules=[], sources=[], gmx='gmx', loaders=[], ignloaderr = True):
    global cmd
    if ':' not in work_dir:
        cmd = local.Local(modules=modules, sources=sources, gmx="gmx", directory=work_dir, loaders=loaders, ignloaderr=ignloaderr)
    else:
        r = work_dir.split(':')
        
        cmd = remote.SSH(r[0], modules=modules, sources=sources, gmx="gmx", directory=r[1], loaders=loaders, ignloaderr=ignloaderr)