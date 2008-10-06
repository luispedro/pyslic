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

import popen2
from numpy.distutils.core import setup, Extension

def readimg_args():
    output,input=popen2.popen2('pkg-config ImageMagick++ --libs')
    input.close()
    tokens=output.readline().split()
    output.close()
    args={ 'libraries'    : [t[2:] for t in tokens if t.startswith('-l')],
           'include_dirs' : [t[2:] for t in tokens if t.startswith('-I')],
           'library_dirs' : [t[2:] for t in tokens if t.startswith('-L')],
    }
    return args

readimg = Extension('readimg', sources = ['readimg.cpp'],  **readimg_args())

if __name__ == '__main__':
    setup (name = 'readimg',
           version = '1.0',
           description = 'Read image using ImageMagick',
           ext_modules = [readimg]
           )
