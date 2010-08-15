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

from __future__ import division, with_statement
import numpy
from scipy.ndimage import gaussian_filter
from ..image import Image, loadedimage

__all__ = [
    'preprocess_collection',
    'FixIllumination',
    'FixIlluminationHVGradient',
    'FixIlluminationRadialGradient',
    'FixIlluminationHVRadialGradient',
    'ConcatPreprocessors',
    'NullPreprocessor',
    ]

def preprocess_collection(imgs,P):
    '''
    P = process_collection(imgs,P)
    
    This function does:

        for img in imgs:
            with loadedimage(img):
                P.see(img)
        P.finish()
        return P
    '''
    for img in imgs:
        with loadedimage(img):
            P.see(img)
    P.finish()
    return P


def _smooth_S(S,grow_i=True,grow_j=True,d_ij=True,d_ij2=True):
    from numpy import c_
    X=numpy.ones_like(S).ravel()
    if grow_i:
        M=numpy.fromfunction(lambda i,j: i, S.shape)
        X=c_[X,M.ravel()]
    if grow_j:
        M=numpy.fromfunction(lambda i,j: j, S.shape)
        X=c_[X,M.ravel()]
    if d_ij or d_ij2:
        ci,cj=S.shape
        ci /= 2
        cj /= 2
        M=numpy.fromfunction(lambda i,j: (i-ci)**2+(j-cj)**2,S.shape)
        if d_ij2:
            X=c_[X,M.ravel()]
        if d_ij:
            M=numpy.sqrt(M)
            X=c_[X,M.ravel()]
    # I would prefer to use something like the Lasso below.
    W=numpy.linalg.lstsq(X,S.ravel())[0]
    S_flat=numpy.dot(X,W)
    S_flat=numpy.reshape(S_flat,S.shape)
    return S_flat


class FixIllumination(object):
    '''
    Fix illumination by computing the average illumination at each pixel in a collection of images
    and then dividing the pixel values by that amount.
    '''
    __slots__ = ['S','channel']
    def __init__(self,channel='protein'):
        self.S=None
        self.channel=channel

    def see(self,img):
        '''
        self.see(img)

        Collect statistics on one image
        '''
        img.lazy_load()
        P=img.channeldata[self.channel]
        if self.S is None:
            self.S = numpy.zeros(P.shape,numpy.float96)
        self.S += P

    def __getstate__(self):
        return (self.S,self.channel)

    def __setstate__(self,state):
        self.S,self.channel = state

    def finish(self):
        '''
        self.finish()

        Signal to the object that all images have been see()n

        @see see
        '''
        assert self.S is not None
        Smin = self.S.min()
        if Smin == 0:
            Smin = 1
        self.S /= Smin
        # float96 is not always very well supported and we no longer need to sum up lots of numbers
        self.S = numpy.array(self.S,float)

    def process(self,img):
        '''
        self.process(img)

        Fix illumination of img
        '''
        img.lazy_load()
        P=img.channeldata[self.channel]
        P /= self.S
        img.channeldata[self.channel] = P

class FixIlluminationRadialGradient(FixIllumination):
    '''
    This is a collection processor that models an illumination
    gradient from the centre of the image outwards.
    '''
    def finish(self):
        FixIllumination.finish(self)
        self.S=_smooth_S(self.S,grow_i=False,grow_j=False,d_ij=True,d_ij2=True)

class FixIlluminationHVGradient(FixIllumination):
    '''
    This is a collection processor that models the illumination
    uneveness as a combination of horizontal and vertical gradient.
    '''
    def finish(self):
        FixIllumination.finish(self)
        self.S=_smooth_S(self.S,grow_i=True,grow_j=True,d_ij=False,d_ij2=False)

class FixIlluminationHVRadialGradient(FixIllumination):
    '''
    This is a collection processor that models the illumination
    uneveness as a combination of horizontal, vertical and radial gradients.
    
    @see FixIlluminationHVGradient 
    @see FixIlluminationRadialGradient
    '''
    def finish(self):
        FixIllumination.finish(self)
        self.S=_smooth_S(self.S,grow_i=True,grow_j=True,d_ij=True,d_ij2=True)

class FixIlluminationGaussianFilter(FixIllumination):
    '''
    This is a collection processor that models the illumination
    uneveness as a combination of horizontal and vertical gradient.
    '''
    __slots__ = ['sigma']
    def __init__(self,sigma=2):
        FixIllumination.__init__(self)
        self.sigma=sigma

    def __getstate__(self):
        base=FixIllumination.__getstate__(self)
        return base+(self.sigma,)

    def __setstate__(self,state):
        FixIllumination.__setstate__(state[:-1])
        self.sigma = state[-1]

    def finish(self):
        FixIllumination.finish(self)
        if self.sigma > 0:
            self.S = gaussian_filter(self.S,self.sigma)
            self.S /= self.S.min()

class ConcatPreprocessors(object):
    __slots__ = ['preprocessors']
    def __init__(self,*preprocessors):
        self.preprocessors = preprocessors

    def see(self,img):
        '''
        self.see(img)

        Calls see() on alls processors in self.preprocessors
        '''
        for P in self.preprocessors:
            P.see(img)

    def finish(self):
        '''
        self.finish(img)

        Calls finish() alls processors in self.preprocessors
        '''
        for P in self.preprocessors:
            P.finish()

    def process(self,img):
        '''
        self.process(img)

        Calls process(img) alls processors in self.preprocessors
        '''
        for P in self.preprocessors:
            P.process(img)

    def __getstate__(self):
        return self.preprocessors

    def __setstate__(self,state):
        self.preprocessors = state

class NullPreprocessor(object):
    def __init__(self):
        pass
    def process(self,img):
        '''
        self.process(img)

        Does nothing
        '''
        pass
# vim: set ts=4 sts=4 sw=4 expandtab smartindent:
