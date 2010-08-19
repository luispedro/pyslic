# -*- coding: utf-8 -*-
# Copyright (C) 2008-2010 Murphy Lab
# vim: set ts=4 sts=4 sw=4 expandtab smartindent:
# Carnegie Mellon University
# 
# Original (Matlab) version by M.V. Boland (10 Aug 98)
# Python version by Luis Pedro Coelho <lpc@cmu.edu>
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

import numpy as np
from mahotas import center_of_mass

__all__ = ['imgmoments','imgcentmoments']
def imgmonents(img, x, y):
    """
     M = imgmonents(IMG, X, Y) calculates the moment MXY for IMAGE
     imgmoments(img, x, y), 
        where IMG is the image to be processed and X and Y define
        the order of the moment to be calculated. For example, 

        imgmoments(image,0,1) calculates the first order moment 
        in the y-direction, and 

        imgmoments(image,0,1)/imgmoments(image,0,0) is the 
        'center of mass (fluorescence)' in the y-direction
    """
    return imgcentmoments(img, x, y, 0, 0)

def imgcentmoments(img,x,y,cofy=None,cofx=None):
    """
    M_xy = imgcentmoments(img,x,y, cofy=None, cofx=None)

    @param cofy and cofx are optional and computed from the image if not given
    """
    if cofy is None or cofx is None:
        print 'calling center_of_mass'
        cofy, cofx = center_of_mass(img)
    if not np.issubdtype(img.dtype, float):
        img = img.astype(float)
    r,c = img.shape
    p = np.arange(c)
    p -= cofx
    p **= x
    inter = np.dot(img, p)
    p = np.arange(r)
    p -= cofy
    p **= y
    return np.dot(inter, p)

