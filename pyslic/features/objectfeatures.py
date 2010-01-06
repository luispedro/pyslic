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
import numpy as np
import warnings
from numpy import *
from scipy import ndimage

from mahotas.bbox import bbox

from ..image import Image
from ..imageprocessing.bweuler import bweuler
from .hullfeatures import hullfeatures
from .mmthin import mmthin
from ..imageprocessing.convexhull import convexhull
from .imgskelfeats import find_branch_points

try:
    import ncreduce
    fast_sum = ncreduce.sum
except ImportError:
    fast_sum = np.sum

def objectfeatures(img):
    '''
    values=objectfeatures(img)

    This implements the object features described in
    "Object Type Recognition for Automated Analysis of Protein Subcellular Location"
    by Ting Zhao, Meel Velliste, Michael V. Boland, and Robert F. Murphy
    in IEEE Transaction on Image Processing
    '''

    protimg = img.channeldata[Image.procprotein_channel]
    dnaimg = img.channeldata.get(Image.procdna_channel,None)
    assert dnaimg is None or protimg.shape == dnaimg.shape, \
        'pymorph.objectfeatures: DNA image is not of same size as Protein image.'

    labeled,N = ndimage.label(protimg,ones((3,3)))

    sofs = np.zeros((N,11))
    if dnaimg is not None:
        dnacofy,dnacofx = ndimage.center_of_mass(dnaimg)
        bindna = (dnaimg > 0)
    
    # According to the documentation, it shouldn't matter if indices is None,
    # but in my version of scipy.ndimage, you *have* to use indices.
    indices = np.arange(1,N+1)
    centers = ndimage.center_of_mass(protimg, labeled, indices)
    for obji in xrange(N):
        binobj = (labeled == (obji+1))
        min1,max1,min2,max2 = bbox(binobj)
        if min1 > 0: min1 -= 1
        if min2 > 0: min2 -= 1
        binobjc = binobj[min1:max1+1,min2:max2+1] # leave a small margin for bweuler()
        protobj = protimg[min1:max1+1,min2:max2+1]
        objimg = protimg * binobj
        cofy,cofx = centers[obji]
        binskel = mmthin(binobjc)
        objhull = convexhull(binobjc)
        no_of_branch_points = fast_sum(find_branch_points(binskel))
        hfeats = hullfeatures(binobjc,objhull)

        sofs[obji,0] = fast_sum(binobjc)
        if dnaimg is not None:
            sofs[obji,1] = np.sqrt((cofy-dnacofy)**2+(cofx-dnacofx)**2)
            sofs[obji,2] = fast_sum(binobj&bindna)/sofs[obji,0]
        sofs[obji, 3] = hfeats[2]
        sofs[obji, 4] = bweuler(binobjc)
        sofs[obji, 5] = hfeats[1]
        sofs[obji, 6] = fast_sum(binskel)
        sofs[obji, 7] = hfeats[0]
        sofs[obji, 8] = sofs[obji,6]/sofs[obji,0]
        sofs[obji, 9] = fast_sum(objimg)/fast_sum(binskel*protobj)
        sofs[obji,10] = no_of_branch_points/sofs[obji,6]
    return sofs

# vim: set ts=4 sts=4 sw=4 expandtab smartindent:
