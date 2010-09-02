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
import warnings
from numpy import *
from scipy import ndimage

from mahotas.bbox import bbox

from ..image import Image
from mahotas import euler
from .hullfeatures import hullfeatures
from mahotas import thin
from mahotas.polygon import fill_convexhull as convexhull
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

    protimg = img.get('procprotein')
    dnaimg = img.channeldata.get('procdna',None)
    assert dnaimg is None or protimg.shape == dnaimg.shape, \
        'pymorph.objectfeatures: DNA image is not of same size as Protein image.'

    labeled,N = ndimage.label(protimg,ones((3,3)))
    if not N:
        return np.zeros((0,11))

    sofs = np.zeros((N,11))
    indices = np.arange(1,N+1)
    if dnaimg is not None:
        dnacofy,dnacofx = ndimage.center_of_mass(dnaimg)
        bindna = (dnaimg > 0)
        # According to the documentation, it shouldn't matter if indices is None,
        # but in my version of scipy.ndimage, you *have* to use indices.
        centers = ndimage.center_of_mass(protimg, labeled, indices)
        if N == 1:
            centers = [centers]
        centers = np.asarray(centers)
        centers -= np.array((dnacofy, dnacofx))
        centers **= 2
        sofs[:,1] = np.sqrt(centers.sum(1))
    locations = ndimage.find_objects(labeled, N)
    sofs[:, 9] = ndimage.measurements.sum(protimg, labeled, indices)
    for obji in xrange(N):
        slice = locations[obji]
        binobj = (labeled[slice] == (obji+1)).copy()
        protobj = protimg[slice]
        binskel = thin(binobj)
        objhull = convexhull(binobj)
        no_of_branch_points = fast_sum(find_branch_points(binskel))
        hfeats = hullfeatures(binobj,objhull)
        sofs[obji,0] = fast_sum(binobj)
        if dnaimg is not None:
            sofs[obji,2] = fast_sum(binobj&bindna[slice])
        sofs[obji, 3] = hfeats[2]
        sofs[obji, 4] = euler(binobj)
        sofs[obji, 5] = hfeats[1]
        sofs[obji, 6] = fast_sum(binskel)
        sofs[obji, 7] = hfeats[0]
        sofs[obji, 9] /= fast_sum(binskel*protobj)
        sofs[obji,10] = no_of_branch_points
    sofs[:,2] /= sofs[:,0]
    sofs[:,8] = sofs[:,6]/sofs[:,0]
    sofs[:,10] /= sofs[:,6]
    return sofs

# vim: set ts=4 sts=4 sw=4 expandtab smartindent:
