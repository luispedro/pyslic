# -*- coding: utf-8 -*-
# Copyright (C) 2008  Murphy Lab
# Carnegie Mellon University
# 
# Written by Luís Pedro Coelho <lpc@cmu.edu>
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

import setuptools
import numpy.distutils.core as numpyutils
from readimg_setup import readimg_args

convexhull = numpyutils.Extension('pyslic/imageprocessing/_convexhull', sources = ['pyslic/imageprocessing/convexhull.cpp'])

readimg = numpyutils.Extension('pyslic.image.io.readimg', sources = ['pyslic/image/io/readimg.cpp'], **readimg_args())
packages=setuptools.find_packages()
packages.remove('tests')

def test_pyversion():
    import sys
    maj,min,_,_,_ = sys.version_info
    if (maj,min) < (2,5):
        print "Your Python interpreter is too old for Pyslic.\nUpgrade to 2.5 or newer.\n"
        sys.exit(1)

test_pyversion()

numpyutils.setup(name='PySLIC',
      version='0.4.4',
      description='Subcellular Location Image Classifier',
      author='Murphy Lab',
      author_email='murphy@mcu.edu',
      url='http://murphylab.cbi.cmu.edu/',
      packages=packages,
      ext_modules = [convexhull,readimg],
      test_suite='nose.collector',
      )


# vim: set ts=4 sts=4 sw=4 expandtab smartindent:
