# -*- coding: utf-8 -*-
# Copyright (C) 2008-2009  Murphy Lab
# vim: set ts=4 sts=4 sw=4 expandtab smartindent:
# Carnegie Mellon University
# 
# Written by Luis Pedro Coelho <lpc@cmu.edu>
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published
# by the Free Software Foundation; either version 2 of the License,
# or (at your option) any later version.
#
# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
# 02110-1301, USA.
#
# For additional information visit http://murphylab.web.cmu.edu or
# send email to murphy@cmu.edu

from try_imports import try_imports
from sys import exit
try:
    import setuptools
except:
    print '''
setuptools not found.

On linux, the package is often called python-setuptools'''
    exit(1)
try_imports()

import numpy.distutils.core as numpyutils

ext_modules = []

packages=setuptools.find_packages()
if 'tests' in packages: packages.remove('tests')

def test_pyversion():
    import sys
    maj,min,_,_,_ = sys.version_info
    if (maj,min) < (2,5):
        print "Your Python interpreter is too old for Pyslic.\nUpgrade to 2.5 or newer.\n"
        sys.exit(1)

test_pyversion()

numpyutils.setup(name='PySLIC',
      version='0.6',
      description='Subcellular Location Image Classifier',
      author='Murphy Lab',
      author_email='murphy@cmu.edu',
      url='http://murphylab.cbi.cmu.edu/',
      packages=packages,
      ext_modules = ext_modules,
      test_suite='nose.collector',
      )


