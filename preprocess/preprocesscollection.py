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
from scipy.ndimage import gaussian_filter
from ..image import Image

__all__ = ['FixIllumination','preprocess_collection']

def preprocess_collection(imgs,P,unload=True):
    '''
    P = process_collection(imgs,P,unload=True)
    
    for img in imgs:
        P.see(img)
        if unload: img.unload()
    P.finish()
    return P

    @param unload: if True, images are unloaded after processing
    '''
    for img in imgs:
        P.see(img)
        if unload:
            img.unload()
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
    __slots__ = ['S','channel','sigma']
    def __init__(self,channel=Image.protein_channel,sigma=2):
        self.S=None
        self.channel=channel
        self.sigma=sigma

    def see(self,img):
        img.lazy_load()
        P=img.channeldata[self.channel]
        if self.S is None:
            self.S = numpy.zeros(P.shape,numpy.float96)
        self.S += P

    def finish(self):
        assert self.S is not None
        self.S /= self.S.min()
        self.S = numpy.array(self.S,float)
        self.S = gaussian_filter(self.S,self.sigma)
        self.S /= self.S.min()

    def process(self,img):
        img.lazy_load()
        P=img.channeldata[Image.protein_channel]
        P /= self.S
        img.channeldata[Image.protein_channel] = P

class FixIlluminationRadial(FixIllumination):
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

# vim: set ts=4 sts=4 sw=4 expandtab smartindent:
