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
from numpy import *
from scipy.ndimage import distance_transform_edt,label, median_filter
from mahotas.thresholding import otsu
import nucleidetection

__all__ = ['voronoi','gvoronoi']
def voronoi(img,centers,distance='euclidean'):
    '''
    labeled = voronoi(img,centers,distance='euclidean')

    distance can be one of
        'euclidean' : use euclidean distance (default)
        'manhatan'  : use manhatan distance
    '''
    labels=zeros_like(img)
    r,c=img.shape
    def seuclidean(y,x,cy,cx):
        dy=y-cy
        dx=x-cx
        return dy*dy+dx*dx
    def manhatan(y,x,cy,cx):
        dy=y-cy
        dx=x-cx
        return abs(dy)+abs(dx)
    dist=seuclidean
    if distance == 'manhatan' or distance == 'cityblock' or distance == 'city_block':
        dist=manhatan
    for y in xrange(r):
        for x in xrange(c):
            best=inf
            for i,(cy,cx) in enumerate(centers):
                now=dist(y,x,cy,cx)
                if now < best:
                    best=now
                    besti=i+1
            labels[y,x]=besti
    return labels

def gvoronoi(labeled,distance='euclidean'):
    """
    segmented = gvoronoi(labeled,distance='euclidean')

    Generalised Voronoi Transform

    INPUT:
    * labeled: an array, of a form similar to the return of scipy.ndimage.label()
    * distance: one of 'euclidean', 'manhatan'

    RETURN
    segmented is of the same size and type as labeled and
        segmented[y,x] is the label of the object at position y,x
    """
    if distance == 'euclidean':
        L1,L2=distance_transform_edt(labeled== 0, return_distances=False,return_indices=True)
    else:
        raise Exception('gvoronoi: Distance "%s" not implemented' % distance)
    return labeled[L1,L2]

def gvoronoi_dna(dnaimg,distance='euclidean'):
    """
    labeled = gvoronoi_dna(dnaimg,distance)

    Generalised Voronoi Transform
    """
    labelednuclei,_=nucleidetection.labelnuclei(dnaimg)
    return gvoronoi(labelednuclei,distance)

# vim: set ts=4 sts=4 sw=4 expandtab smartindent:
