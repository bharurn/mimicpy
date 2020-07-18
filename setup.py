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

import mimicpy

setup(
    name='mimicpy',
    version=mimicpy.__version__,
    zip_safe=False,
    description='Python tools to prepare MiMiC QM/MM runs.',
    author=mimicpy.__authors__,
    author_email='b.raghavan@fz-juelich.de',
    packages=get_package(),
    #install_requires=get_reqs(),
    entry_points = {
        'console_scripts': ['mimicpy = mimicpy.__main__:main'],
    })