# This is part of MiMiCPy

"""

This module contains functions to run commands on local shell

"""

import subprocess
import os
from ..utils.errors import defaultHook

decoder = 'utf-8'

def run(cmd, shell_ex, remove_from_out, stdin=None, hook=None):
    if stdin == None or stdin.strip() == '':
        p = subprocess.Popen(cmd, shell=True, executable=shell_ex, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        out = p.stdout.read().decode(decoder)
    else:
        p = subprocess.Popen(cmd, shell=True, executable=shell_ex, stdin=subprocess.PIPE,\
                             stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        
        inp = stdin
        if isinstance(stdin, str):
            inp = stdin.encode(decoder)
        
        out = p.communicate(input=inp)[0].decode()
    
    out = out.replace(remove_from_out, '')
    
    if not hook: hook = defaultHook
    
    if hook: hook(cmd.split(';')[-1], out)
    
    return out

def runbg(self, cmd, shell_ex, remove_from_out, hook=None):
    subprocess.Popen(cmd, shell=True, executable=shell_ex, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    # add communication till output not changes like remote

def ls(dirc, file_eval, dir_eval):
        
    files = []
        
    if dirc:
        dirs = os.listdir(dirc)
    else:
        dirs = os.listdir()
        
    for file in dirs:
        if os.path.isdir(file):
            if dir_eval(file): files.append(file)
        elif file_eval(file):
            files.append(file)
            
    return files
