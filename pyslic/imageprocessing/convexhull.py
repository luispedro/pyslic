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
from ..imageprocessing.bbox import bbox
from numpy import *
import Image, ImageDraw
from scipy.misc import fromimage
import warnings

__all__ = ['convexhull','computeConvexHull']

def convexhull(bwimg):
    """
    hull = convexhull(bwimg)

    Compute the convex hull of the binary image and return it as a binary mask
    """
    Y,X=where(bwimg)
    P=list(zip(Y,X))
    if len(P) <= 3:
        return bwimg
    H=computeConvexHull(P)
    
# This is kind of sucky code, but it works:
    im=Image.new('L',bwimg.shape,0)
    draw=ImageDraw.Draw(im)
    draw.fill=1
    draw.polygon(H,1,1)
    return fromimage(im).T

def _isLeft(p0,p1,p2):
    """
    tests if a point is Left|On|Right of an infinite line.
    Input:  three points p0, p1, and p2

    Return:
        >0 for p2 left of the line through p0 and p1
        =0 for p2 on the line
        <0 for p2 right of the line
    """
    # Copyright 2001, softSurfer (www.softsurfer.com)
    # This code may be freely used and modified for any purpose
    # providing that this copyright notice is included with it.
    # SoftSurfer makes no warranty for this code, and cannot be held
    # liable for any real or imagined damage resulting from its use.
    # Users of this code must verify correctness for their application.
    return (p1[0]-p0[0])*(p2[1]-p0[1]) - (p2[0]-p0[0])*(p1[1]-p0[1])

def _inPlaceScan(P,reverse):
    P.sort(reverse=reverse)
    h=1
    N=len(P)
    for i in xrange(1,N):
        while h >= 2 and _isLeft(P[h-2],P[h-1],P[i]) >= 0:
            h -= 1
        P[i],P[h]=P[h],P[i]
        h += 1
    return h
def computeConvexHull(P):
    """
    From ``Space-Efficient Planar Convex Hull Algorithms''
        by Bronnimann et al.
    """
    try:
        import _convexhull
        return _convexhull.computeConvexHull(P)
    except ImportError:
        warnings.warn('C code for convex hull failed. Resorting to (slow) python (Error: %s)' % e)
        h=_inPlaceScan(P,False)
        for i in xrange(h-1):
            P[i],P[i+1]=P[i+1],P[i]
        P2=P[h-2:]
        h_=_inPlaceScan(P2,True)
        return P[:h-2]+P2[:h_]

# vim: set ts=4 sts=4 sw=4 expandtab smartindent:
