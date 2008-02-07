## This file was originally part of the octave-forge project
## Ported to python by Luis Pedro Coelho <luis@luispedro.org> (February 2008)
## Copyright (C) 2006  Soren Hauberg
## Copyright (C) 2008  Luis Pedro Coelho (Python port)
## 
## This program is free software; you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation; either version 2, or (at your option)
## any later version.
## 
## This program is distributed in the hope that it will be useful, but
## WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
## General Public License for more details. 
## 
## You should have received a copy of the GNU General Public License
## along with this file.  If not, see <http://www.gnu.org/licenses/>.

from numpy import *

def bwperim(bw,n=4):
    """
    Find the perimeter of objects in binary images.

    A pixel is part of an object perimeter if its value is one and there
    is at least one zero-valued pixel in its neighborhood.

    By default the neighborhood of a pixel is 4 nearest pixels, but
    if @var{n} is set to 8 the 8 nearest pixels will be considered.
    """

    if n != 4 and n != 8:
        raise Exception('bwperim: n must be 4 or 8')
    out = bw.copy()
    rows,cols = bw.shape

    # Translate image by one pixel in all directions
    north = zeros((rows,cols))
    south = zeros((rows,cols))
    west = zeros((rows,cols))
    east = zeros((rows,cols))

    north[:-1,:] = bw[1:,:]
    south[1:,:]  = bw[:-1,:]
    west[:,:-1]  = bw[:,1:]
    east[:,1:]   = bw[:,:-1]
    if n == 4:
        idx = (north == bw) & (south == bw) & (west == bw) & (east == bw);
    else:
        north_east = zeros((rows, cols))
        north_west = zeros((rows, cols))
        south_east = zeros((rows, cols))
        south_west = zeros((rows, cols))
        north_east[0:-1, 1:]   = bw[1:, :-1]
        north_west[0:-1, 1:-1] = bw[1:, 1:]
        south_east[1:, 1:]     = bw[:-1, :-1]
        south_west[1:, 1:-1]   = bw[:-1, 1:]
        idx = (north == bw) & (north_east == bw) & \
              (east  == bw) & (south_east == bw) & \
              (south == bw) & (south_west == bw) & \
              (west  == bw) & (north_west == bw);
    out[idx] = 0;
    return out


# vim: set ts=4 sts=4 sw=4 expandtab smartindent:
