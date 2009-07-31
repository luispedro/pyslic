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
    angle = (2*np.pi)/points
    feature = []
    lis = []
    def binary(i, n):
        j = 0
        return list((0,1)[i>>j & 1] for j in xrange(n-1,-1,-1))
    for q in xrange(2**points):
        lis = binary(q,q)
        if len(lis) < points:
            for i in xrange(points-len(lis)):
                lis.insert(0,0)
        if len(lis) > points:
            for i in xrange(len(lis)-points):
                del lis[0]
        feature.append(lis)
        lis = []
    final = np.array([0] * (len(feature)))
    
    for row in xrange(radius, image.shape[0]-radius-1):
        for pix in xrange(radius, image.shape[1]-radius-1):
            sign = [] # For the comparison function
            center = image[row,pix]
     
            for l in xrange(points): # For every point calculate these features
                x = radius * np.sin(angle*l)
                y = radius * np.cos(angle*l)
                ix = int(x)
                iy = int(y)
                a = image[(row + iy),(pix + ix)]
                b = image[(row + iy + 1),(pix + ix)]
                c = image[(row + iy),(pix + ix +1)]
                d = image[(row + iy+1),(pix + ix+1)]
                dx = x-ix
                dy = y-iy
                e = a+(c-a)*(dy)
                f = b+(d-b)*(dy)
                r = e+(f-e)*(dx)
                sign.append(int(center > r))
            sign = np.array(sign)
            bestval = sign.copy()

            for n in xrange(len(sign)):
                cur = np.roll(sign, n)
                for curbit, bestbit in zip(cur,bestval):
                    if curbit != bestbit:
                        if curbit < bestbit:
                            bestval = cur
                        break
            bestval = list(bestval)
            g = feature.index((bestval))
            final[g] += 1
    return final

