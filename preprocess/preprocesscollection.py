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

    def process(self,img):
        img.lazy_load()
        P=img.channeldata[Image.protein_channel]
        P /= self.S
        img.channeldata[Image.protein_channel] = P

# vim: set ts=4 sts=4 sw=4 expandtab smartindent:
