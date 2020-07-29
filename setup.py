#!/usr/bin/env python

import os
import sys
from setuptools import setup, find_packages

def get_package():
    root = 'mimicpy'
    return  [root]+[root+'.'+i for i in find_packages(root)]

def get_reqs():
    with open('requirements.txt', 'r') as f:
        reqs = f.read().splitlines()
    
    return reqs

def get_details(detail, deflt):
    path = 'mimicpy/_'+detail
    if not os.path.isfile(path):
        return deflt
        
    with open(path) as f:
        txt = f.read()
        if len(txt.splitlines()) == 1 and '=' in txt:
            return f.read().split('=')[1].strip()
        else:
            return deflt

setup(
    name='mimicpy',
    version=get_details('version', 1.0),
    zip_safe=False,
    description='Python tools to prepare MiMiC QM/MM runs.',
    author=get_details('authors', "--"),
    author_email='b.raghavan@fz-juelich.de',
    packages=get_package(),
    install_requires=get_reqs(),
    entry_points = {
        'console_scripts': ['mimicpy = mimicpy.__main__:main'],
    })
