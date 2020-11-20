#!/usr/bin/env python

import sys
import os
import platform
from shutil import copyfileobj
from setuptools import setup, find_packages
from setuptools.command.install import install
from setuptools.command.develop import develop
from setuptools.command.egg_info import egg_info

def get_package():
    root = 'mimicpy'
    return  [root]+[root+'.'+i for i in find_packages(root)]

class PostBaseCommand(object):
    """Base class for post-installation code"""

    def copy_script(self, source, dest):
        if not os.path.isfile(source):
            print("Cannot find source file {}".format(source))
            sys.exit(1)
        else:
            with open(source, 'r') as fsrc, open(dest, 'a') as fdest:
                copyfileobj(fsrc, fdest)

    def run(self):
        super().run()
        pymol_dir = os.environ.get('PYMOLDIR', default=None)
        vmd_dir = os.environ.get('VMDDIR', default=None)

        if pymol_dir:
            pymolrc = os.path.join(pymol_dir, '.pymolrc.py')
            self.copy_script('plugins/pymol.py', pymolrc)
            print("Wrote PyMOL settings to file {}".format(pymolrc))

        if vmd_dir:
            if platform.system() == 'Windows':
                vmdrc = os.path.join(vmd_dir, 'vmd.rc')
            else:
                vmdrc = os.path.join(vmd_dir, '.vmdrc')

            self.copy_script('plugins/vmd.tcl', vmdrc)
            print("Wrote VMD settings to file {}".format(vmdrc))

class PostInstallCommand(PostBaseCommand, install):
    """Post-installation code for installation mode"""
    pass

class PostDevelopCommand(PostBaseCommand, develop):
    """Post-installation code for develop mode"""
    pass

class PostEggInfoCommand(PostBaseCommand, egg_info):
    """Post-installation code for egg info mode"""
    pass

with open("README.md", "r") as f:
    long_description = f.read()

setup(
    name='mimicpy',
    version='1.0',
    zip_safe=True,
    description='Python tools to prepare MiMiC QM/MM runs.',
    author="Bharath Raghavan and Florian Schackert",
    author_email='b.raghavan@fz-juelich.de',
    long_description=long_description,
    long_description_content_type="text/markdown",
    license="GNU General Public License v3.0",
    platforms="OS Independent",
    url="https://github.com/bharurn/mimicpy",
    packages=get_package(),
    install_requires=['numpy>=1.12.0', 'pandas>=0.24.0'],
    python_requires='>=3.5',
    classifiers=[
            "Programming Language :: Python :: 3",
            "License :: OSI Approved :: GNU License",
            "Operating System :: OS Independent",
        ],
    entry_points = {
        'console_scripts': [
            'mimicpy = mimicpy.__main__:main',
            'mimicpy_vmd = mimicpy.__main_vmd__:main'
        ],
    },
    cmdclass={
            'install': PostInstallCommand,
            'develop': PostDevelopCommand,
            'egg_info': PostEggInfoCommand,
    })
