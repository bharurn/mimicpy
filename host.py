from pygmx.shell import remote, local

cmd = local.Local()
gmx = 'gmx'

def setup(work_dir, modules=[], sources=[], loaders=[], ignloaderr = True):
    global cmd
    work_dir = work_dir.strip()
    
    if isinstance(sources, str): sources = [sources]
    if isinstance(modules, str): modules = [modules]
    if isinstance(loaders, str): loaders = [loaders]
    
    if ':' not in work_dir:
        print("Setting local machine as host")
        
        cmd = local.Local(modules=modules, sources=sources, directory=work_dir, loaders=loaders, ignloaderr=ignloaderr)
    else:
        
        r = work_dir.split(':')
        print(f"Setting remote machine {r[0]} as host..")
        
        cmd = remote.SSH(r[0], modules=modules, sources=sources, directory=r[1], loaders=loaders, ignloaderr=ignloaderr)
        
def close():
    if isinstance(cmd, remote.SSH):
        cmd.__del__()
        print(f"Closed connections to {cmd.name}..")
    else:
        raise Exception("Cannot close connection! Currently running on localhost.")