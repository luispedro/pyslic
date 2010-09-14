# -*- coding: utf-8 -*-
# Copyright (C) 2008  Murphy Lab
# Carnegie Mellon University
# vim: set ts=4 sts=4 sw=4 expandtab smartindent:
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
from mahotas.morph import majority_filter
from ..imageprocessing.thresholding import threshold
from ..image import Image, loadedimage
from scipy import ndimage

def threshold_segment(dna, threshold_method='otsu', smooth=None, median_size=5, min_obj_size=2500):
    '''
    labeled = threshold_method(dna, threshold_method='otsu',median_size=5,min_obj_size=2500)

    Simple threshold-based segmentation

    Params
    ------
        * dna: either a pyslic.Image or a DNA image
        * threshold_method: thresholding method to use. Can be either a function or 
            string which denotes the name of a function in pyslic.imageprocessing.thresholding (default: otsu)
        * smooth: either None (no smoothing) or a sigma value for a gaussian blur (default: None)
        * median_size: median filter size (default: 5). Set to None to skip filtering.
        * min_obj_size: minimum object size (default: 2500)
    '''
    if type(dna) == Image:
        with loadedimage(dna):
            return threshold_segment(dna.get('dna'),threshold_method,smooth, median_size,min_obj_size)
    if smooth is not None:
        dna = ndimage.gaussian_filter(dna,smooth)
    T = threshold(dna,threshold_method)
    binimg = dna > T
    if median_size is not None:
        binimg = majority_filter(binimg, median_size)
    L,N = ndimage.label(binimg)
    if N == 0:
        return L
    sizes = numpy.array(ndimage.sum(binimg,L,numpy.arange(N+1)))
    for oid in numpy.where(sizes < min_obj_size)[0]:
        L[L == oid] = 0
    L,N = ndimage.label(L != 0)
    return L


