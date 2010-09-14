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
from scipy.ndimage import histogram
from mahotas.thresholding import rc, otsu

__all__=['threshold', 'rc','murphy_rc','otsu','softthreshold','hardthreshold']

def threshold(img,thresh):
    '''
    T = threshold(img, thresh)

    thresh can be:
        * None: returns -1
        * a number: returns thresh
        * a function: returns thresh(img)
        * a string:
            one of ('otsu','rc','murphy_rc','mean')
    '''
    if thresh is None:
        return -1
    if type(thresh) is str:
        if thresh == 'otsu':
            return otsu(img)
        if thresh == 'rc':
            return rc(img)
        if thresh == 'murphy_rc':
            return murphy_rc(img)
        if thresh == 'mean':
            return img.mean()
        raise ValueError("pyslic.threshold: Cannot handle argument '%s'" % thresh)
    if callable(thresh):
        return thresh(img)
    return thresh

def softthreshold(img,T):
    '''
    softthreshold(img,T)


    Implement a soft threshold:
        img[i] = max(img[i]-T,0)

    Processes the image inplace, return a reference to img.

    Use
    B = softthreshold(A.copy(),T)
    to get a copy.

    @see hardthreshold
    '''
    img -= np.minimum(img,T)
    return img

def hardthreshold(img,T):
    '''
    hardthreshold(img,T)


    Implement a soft threshold:
        img[i] = (img[i] if img[i] > T else 0)

    Processes the image inplace, return a reference to img.

    Use
    B = hardthreshold(A.copy(),T)
    to get a copy.

    @see softthreshold
    '''
    img *= (img > T)
    return img


def murphy_rc(img,ignore_zeros=False):
    """
    T = murphy_rc(img)
    
    Calculate a threshold according to Murphy's adaptation of the RC method.

    @param ignore_zeros: Whether to ignore zero valued pixels (default: False)
        Murphy's Matlab implementation always ignores zero valued pixels.
    """
    pmax = img.max()
    return pmax-rc(pmax-img, ignore_zeros=ignore_zeros)


# vim: set ts=4 sts=4 sw=4 expandtab smartindent:
