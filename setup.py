#!/usr/bin/env python
# -*- coding: ISO-8859-1 -*-
# $Id: setup.py 77 2009-05-30 18:22:58Z cristian_docking $
#
#  Copyright (C) 2009   Cristian S. Rocha (crocha@dc.uba.ar)
#  Licensed to PSF under a Contributor Agreement.
"""Distutils based setup script for Molecular Toolkit.

This uses Distutils (http://python.org/sigs/distutils-sig/) the standard
python mechanism for installing packages. For the easiest installation
just type the command (you'll probably need root privileges for that):

    python setup.py install

This will install the library in the default location. For instructions on
how to customize the install procedure read the output of:

    python setup.py --help install

In addition, there are some other commands:

    python setup.py clean -> will clean all trash (*.pyc and stuff)
    python setup.py test  -> will run the complete test suite
    python setup.py test_core -> will run only tests concerning core features
    python setup.py test_doc -> will run tests on the examples of the documentation

To get a full list of avaiable commands, read the output of:

    python setup.py --help-commands

Or, if all else fails, feel free to write to the sympy list at
mtk@googlegroups.com and ask for help.
"""

#from setuptools import setup, find_packages, Extension
import setuptools
from numpy.distutils.core import Extension
#from numpy.distutils.core import setup
from setuptools import find_packages, setup
from mtk.__version__ import version
import sys
import mtk

from mtk import have_gpu, have_mayavi, have_vtk

# Check arquitecture
debug=1
if debug:
    ext_c_conf={
        'extra_compile_args':['-g3','-gdwarf-2'],
    }
else:
    ext_c_conf={
        'extra_compile_args': ['-mtune=native', '-O3','-DNDEBUG'],
        'define_macros':[('NDEBUG', '1'),],
    }

# Make sure I have the right Python version.
if sys.version_info[1] < 4:
    print "MTK requires Python 2.4 or newer. Python %d.%d detected" % \
          sys.version_info[:2]
    sys.exit(-1)


# List of requires
requires = [
#    'ply',
    'numpy>=1.1.0',
]

# Check some alternatives requirements
if have_gpu():
    requires.append('pycuda>=0.94rc'),
    requires.append('jinja2>=2.5')
else:
    print "Enable Nvidia CUDA code installing pycuda and jinja2."

if have_mayavi():
    requires.append('Mayavi>=3.3.0')
else:
    print "Enable Visualization installing mayavi2"

setup(
    name='mtk',
    version='.'.join(map(str, version)),
    packages=find_packages(),
    scripts=['scripts/drawmol',
        'scripts/drawvol',
        'scripts/drawcurvature',
        'scripts/genvol',
        'scripts/prm2csv',
        'scripts/rotate',
        'scripts/dock',
        'scripts/mtktool',
        'scripts/dx2npy',
        'scripts/cazalsconnolling',
        'scripts/fftconnolling',
        'scripts/crocha',

        'scripts/ply2vtp',
        'scripts/vtp2vti',
        'scripts/fillvti',
        'scripts/vti2vtp',
        'scripts/runexp',

#       'scripts/intersection',
#       'scripts/volhist',
#       'scripts/volrmsd',
            ],

    install_requires = requires,

    package_data = {
        'mtk': ['data/csv/*.csv', 'data/pdb/*.pdb'],
        },
    ext_modules = [ Extension('mtk.geometry.vol_c',
                              ['src/vol_c.c'],
                              depends=['src/vol_c.h'],
                              **ext_c_conf
                             ),
                    Extension('mtk.energy.electrostatic_c',
                              ['src/electrostatic_c.c'],
                              depends=['src/vol_c.h'],
                              **ext_c_conf
                             ),
                  ],

    # Metadata
    author='Cristian S. Rocha',
    author_email='crocha@dc.uba.ar',
    description='Molecular Toolkit for experimental functions',
    license = "PSF",
    keywords = 'sructural bioinformatics toolkit',
    url='http://www.python.org/',
    # TODO: long_description=''
    # TODO: download_url=''
    # TODO: classifiers=''

    # Test Suite
    test_suite = "mtk.tests.test_all.getTestSuite",
    # Documentation
    #docs_source = "docs",
)

