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
from scipy import ndimage
from mahotas import bwperim
from mahotas.bbox import croptobbox
from mahotas.polygon import fill_convexhull as convexhull
import warnings

__all__ = ['dna_size_shape','border_regions','filter_labeled']

def dna_size_shape(labeled, scale=1., minarea=float('-inf'), maxarea=float('+inf'), minroundness=-1.):
    '''
    positives = dna_size_shape(dnamasks, scale=1., minarea=-Inf, maxarea=+Inf, minroundness=0)

    Only accepts DNA objects that fulfill the following criterion:
            * are greater than minarea
            * are smaller than maxarea
            * are rounder than minroundness
    '''
    nr_objects = labeled .max()
    positives = np.zeros(nr_objects+1, np.bool)
    minarea /= scale
    maxarea /= scale
    for obj in xrange(1,nr_objects+1):
        objimg = croptobbox(labeled == obj)
        area = objimg.sum()
        if area > maxarea or area < minarea:
            continue
        hull = convexhull(objimg)
        hullArea = hull.sum()
        hullPerim = bwperim(hull).sum()
        roundness = hullPerim**2/(4*np.pi*hullArea)
        if roundness < minroundness:
            continue
        positives[obj] = 1
    return positives

def border_regions(mask_or_labeled):
    '''
    positives = border_regions(mask_or_labeled)

    Removes anything that touches the border.
    '''
    nr_objs = mask_or_labeled.max()
    if nr_objs == 1:
        labeled,nr_objs = ndimage.label(mask_or_labeled)
    else:
        labeled = mask_or_labeled
    positives = np.ones(nr_objs+1, np.bool)

    negs = np.unique(labeled[0])
    positives[negs] = False

    negs = np.unique(labeled[-1])
    positives[negs] = False

    negs = np.unique(labeled[:,0])
    positives[negs] = False

    negs = np.unique(labeled[:,-1])
    positives[negs] = False
    return positives

def filter_labeled(labeled, positives):
    '''
    labeled, N = filter_labeled(labeled, positives)

    Performs a filtering version of label(), where object OBJ is kept
        only if positives[OBJ]
    '''
    new_label = np.cumsum(positives)
    new_label[~positives] = 0
    assert labeled.max() < len(new_label), 'pyslic.segmentation.filter_labeled: Positives is too small!'
    labeled = new_label[labeled]
    return labeled,int(positives.sum())

# vim: set ts=4 sts=4 sw=4 expandtab smartindent:
