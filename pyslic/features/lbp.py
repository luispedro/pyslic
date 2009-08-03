# -*- coding: utf-8 -*-
# Copyright (C) 2008-2009  Murphy Lab, Carnegie Mellon University
#
# Written by Robert Webb and Luís Pedro Coelho <lpc@cmu.edu>
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
from scipy import ndimage

from scipy import weave
from scipy.weave import converters
def lbp(image, radius, points):
    '''
    features = lbp(image, radius, points)

    Compute Linear Binary Patterns

    Parameters
    ----------
        * image: is the raw image to be evaluated.
        * radius: radius (in pixels)
        * points: nr of points to consider

    Output
    ------
        * features: histogram of features.
    

    Reference
    ---------

        Gray Scale and Rotation Invariant Texture Classification with Local Binary Patterns
            Ojala, T. Pietikainen, M. Maenpaa, T. LECTURE NOTES IN COMPUTER SCIENCE (Springer)
            2000, ISSU 1842, pages 404-420  
    '''
    image = image.astype(np.float)
    angles = np.arange(points) * (2*np.pi)/float(points)
    delta_xs = radius * np.sin(angles)
    delta_ys = radius * np.cos(angles)
    final = np.zeros(2**points)

    def compute_canonical(cur):
        bestval = cur
        for n in xrange(points):
            is_left_bit = (cur & 1)
            cur >>= 1
            if is_left_bit:
                cur |= (1 << points)
            if cur < bestval: bestval = cur
        return bestval

    if points < 20:
        canonical_cache = np.zeros(2**points, np.int32) - 1
        def canonical(input):
            if canonical_cache[input] == -1:
                canonical_cache[input] = compute_canonical(input)
            return canonical_cache[input]
    else:
        canonical = compute_canonical

    coordinates = np.empty( (2,points) )
    rs = np.empty(points, np.float)
    for row in xrange(radius, image.shape[0]-radius):
        for col in xrange(radius, image.shape[1]-radius):
            center = image[row, col]
            coordinates[0] = row
            coordinates[0] += delta_xs
            coordinates[1] = col
            coordinates[1] += delta_ys
            ndimage.interpolation.map_coordinates(image, coordinates, order=1, output=rs)
            code = (2**np.arange(points) * (center > rs)).sum()
            code = canonical(code)
            final[code] += 1
    return final

