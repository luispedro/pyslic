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
import numpy
from scipy import ndimage
from mahotas.morph import majority_filter
from .thresholding import otsu, rc, murphy_rc

def nonzeromin(img):
    '''
    Returns the minimum non zero element in img.
    '''
    A = A.copy().ravel()
    A.sort()
    if A[-1] == 0: return 0
    pos = np.searchsorted(A, 1)
    return A[pos]


def mean_filter(img,size=3):
    '''
    meanimg = mean_filter(img, size=3)

    computes

    meanimg[i,j] = img[i-size//2:i+size//2+1,j-size//2:j+size//2+1].mean()

    i.e., meanimg[i,j] is the mean of the squared centred around (i,j)
    '''
    if len(img.shape) == 2:
        mask=numpy.ones((size,size))/size/size
    elif len(img.shape) == 3:
        mask=numpy.ones((size,size,size))/size/size/size
    else:
        raise ValueError,'mean_filter: img is of wrong shape (can only handle 2d & 3d)'
    return ndimage.convolve(img,mask)

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
        func=ndimage.median_filter
    else:
        raise ArgumentErrorType,"localthresholding: unknown method '%s'" % method
    return img > func(img,size)

def localglobal(img,ignore_zeros=True,globalmethod='otsu',localmethod='mean',localsize=8):
    '''
    Perform both local and global thresholding.
    
    result[i,j] = (img[i,j] > global_threshold) * (img[i,j] > local_threshold[i,j])

    @param img: The image
    @param ignore_zeros: Whether to ignore zero-valued pixels
    @param globalmethod: Global method to use ('otsu', 'rc', or 'murphy_rc')
    @param localmethod: Which local method to use (@see localthresholding)
    @param localsize: Size parameter for local thresholding (@see localthresholding)
    '''
    localobjects=localthresholding(img,method=localmethod,size=localsize)
    if globalmethod == 'otsu':
        T=otsu(img,ignore_zeros=ignore_zeros)
    elif globalmethod == 'rc':
        T=rc(img,ignore_zeros=ignore_zeros)
    elif globalmethod == 'murphy_rc':
        T=murphy_rc(img,ignore_zeros=ignore_zeros)
    else:
        raise ValueError,"localglobal: globalmethod '%s' not recognised." % globalmethod
    globalobjects=img > T
    return localobjects * globalobjects


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
        pmin=nonzeromin(img)
    else:
        pmin=img.min()
    thresholds=pmin+firstThreshold+(img.max()-pmin-firstThreshold)//nrThresholds*numpy.arange(nrThresholds)
    Ts=[majority_filter(img>T) for T in thresholds]
    obj_count = 0
    Ts.append(Ts[0]*0)
    labeled0,N0=ndimage.label(Ts[0])
    for T1 in Ts[1:]:
        labeled1,N1=ndimage.label(T1)
        for obj in xrange(N0):
            binimg=(labeled0 == (obj+1))
            objects1=(labeled1*binimg)
            if not objects1.any() or ndimage.label(objects1)[1] == 1:
                obj_count += 1
                output[binimg]=obj_count
        labeled0=labeled1
        N0=N1

    return output,obj_count

# vim: set ts=4 sts=4 sw=4 expandtab smartindent:
