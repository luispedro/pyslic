# -*- coding: utf-8 -*-
# Copyright (C) 2008  Murphy Lab
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

from __future__ import division
import numpy as np
import numpy
from scipy import ndimage
from numpy import array, asarray, ones, uint32, uint8, zeros
from nhistogram import nhistogram

__all__ = ['fullhistogram','majority_filter','nonzeromin']

def fullhistogram(img):
    """
    H = fullhistogram(img)

    Return a histogram with bins
        0, 1, ..., img.max()
    """
    try:
        from scipy import weave
        histogram = np.zeros(img.max()+1,np.int32)
        img = np.ascontiguousarray(img.ravel())
        N = img.size
        code = '''
#line 46 "basics.py"
        for (int i = 0; i != N; ++i) {
            ++histogram[static_cast<long int>(img[i])];
        }
        '''
        weave.inline(code, ['img','N','histogram'])
        return histogram
    except ImportError:
        maxt = img.max()
        if maxt == 0:
            return np.array([img.size])
        return nhistogram(img, numpy.arange(maxt+2))[0]

def majority_filter(bwimg, N = 3):
    """
    binimg = majority_filter(bwimg, N = 3)

    Implement a majority filter of size N.

    bwimg is interpreted as a binary image (non-zero pixels correspond to one/True).
    Returns a binary image.

    @param N : size of filter (must be an odd integer
    """
    if not N & 1:
        import warnings
        warnings.warn('majority_filter: size argument must be odd.')
        N += 1
    try:
        from scipy import weave
        from scipy.weave import converters
        N=int(N)
        r,c=bwimg.shape
        bwimg=asarray(bwimg,uint8)
        output=zeros((r,c),uint8)
        T=round(N**2/2)
        code = '''
#line 85 "basics.py"
        for (int y = 0; y != r-N+1; ++y) {
            for (int x = 0; x != c-N+1; ++x) {
                int count = 0;
                for (int dy = 0; dy != N; ++dy) {
                    for (int dx = 0; dx != N; ++dx) {
                        if (bwimg(y+dy,x+dx)) ++count;
                    }
                }
                if (count >= T) {
                    output(y+int(N/2),x+int(N/2)) = 1;
                }
            }
        }
        '''
        weave.inline(code,
            ['r','c','bwimg','N','T','output'],
            type_converters=converters.blitz
            )
        return output
    except Exception, e:
        import warnings
        warnings.warn('scipy.weave failed (Error: %s). Resorting to (slow) Python code.' % e)
        assert N < 2 ** 16
        bwimg = asarray(bwimg > 0,uint32)
        filter = ones((N,N))
        bwimg = ndimage.convolve(bwimg,filter)
        return bwimg > (N**2/2)


# vim: set ts=4 sts=4 sw=4 expandtab smartindent:
