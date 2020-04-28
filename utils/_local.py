import subprocess
from os import path
import os

decoder = 'utf-8'

def run(cmd, stdin=None, errorHandle=None):
    args = cmd.split()
    if stdin == None: 
        cmd = subprocess.run(args, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    elif stdin.strip() == '':
        cmd = subprocess.run(args, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    else:
        if isinstance(stdin, bytes):
            cmd = subprocess.run(args, stdout=subprocess.PIPE, stderr=subprocess.PIPE, input=stdin)
        else:
            cmd = subprocess.run(args, stdout=subprocess.PIPE, stderr=subprocess.PIPE, input=stdin.encode('utf-8'))
            
    if errorHandle != None:
        errorHandle(cmd.stdout.decode(decoder))
        errorHandle(cmd.stderr.decode(decoder))
        
    if cmd.returncode != 0:
        raise Exception(f"Error executing {args[0]}!\nCommand executed: {subprocess.list2cmdline(args)}\n{cmd.stdout.decode(decoder)}\n{cmd.stderr.decode(decoder)}")
    
    return cmd.stdout.decode(decoder)
    
def runinSeq(*args, stdin=''):
    val = run(args[0], stdin=stdin)
    
    for i in args[1:-1]:
        if i[0] == 'print':
            print(i[1])
        else:
            val = run(i, stdin=val)
    
    return run(*args[-1], stdin=val)

def fileExists(file):
    if path.isfile(file) or path.isdir(file):
        return True
    else:
        return False

def cd(directory, mkdir=False):
    if directory.strip() == '.':
        return -1
        
    if not fileExists(directory):
        if mkdir:       
            run('mkdir', directory)
            return 1
        else:
            raise Exception(f'Directory {directory} not found')
    else:
        os.chdir(directory)
        return 0