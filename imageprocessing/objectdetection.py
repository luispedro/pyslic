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

from __future__ import division
import numpy
from scipy.ndimage import label, convolve
from .basics import majority_filter


def mean_filter(img,size=3):
    '''
    meanimg = mean_filter(img, size=3)

    computes

    meanimg[i,j] = img[i-size//2:i+size//2+1,j-size//2:j+size//2+1].mean()
    '''
    mask=numpy.ones((size,size))/size/size
    return convolve(img,mask)

def localthresholding(img,method='mean',size=8):
    '''
    thresholded = localthresholding(img,method='mean',size=8)

    @param method: One of
        'mean': compute the mean pixel value
        'median': compute the median pixel value

    @param size size of the window to use
    '''
    if method == 'mean':
        func=mean_filter
    elif method == 'median':
        func=median_filter
    else:
        raise ArgumentErrorType,"localthresholding: unknown method '%s'" % method
    return img > func(img,size)


def multithreshold(img,ignore_zeros=True,firstThreshold=20,nrThresholds=5):
    '''
    labeled,N = multithreshold(img, ignore_zeros = True)

    Performs multi thresholding (which is a form of oversegmentation).

    labeled is of the same size and type as img and contains different labels for 
    the N detected objects (the return of this function is  similar to that of scipy.ndimage.label())

    @param img: The input image
    @param ignore_zeros: Don't take zero pixels into account
    '''
    output=numpy.zeros_like(img)
    if ignore_zeros:
        pmin=img.max()
        r,c=img.shape
        # This is the kind of loop that looks ripe for scipy.weave
        for i in xrange(r):
            for j in xrange(c):
                if img[i,j] and img[i,j] < pmin:
                    pmin = img[i,j]
    else:
        pmin=img.min()
    thresholds=pmin+firstThreshold+(img.max()-pmin-firstThreshold)//nrThresholds*numpy.arange(nrThresholds)
    Ts=[majority_filter(img>T) for T in thresholds]
    obj_count = 0
    Ts.append(Ts[0]*0)
    labeled0,N0=label(Ts[0])
    for T1 in Ts[1:]:
        labeled1,N1=label(T1)
        for obj in xrange(N0):
            binimg=(labeled0 == (obj+1))
            objects1=(labeled1*binimg)
            if not objects1.any() or label(objects1)[1] == 1:
                obj_count += 1
                output[binimg]=obj_count
        labeled0=labeled1
        N0=N1

    return output,obj_count

# vim: set ts=4 sts=4 sw=4 expandtab smartindent:
