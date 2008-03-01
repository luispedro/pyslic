## The original version is part of the octave-forge project (original files: bweuler.m and applylut.m)
## Copyright (C) 2004 Josep Mones i Teixidor
## Copyright (C) 2008 Luis Pedro Coelho - Ported to Python
##
## This program is free software; you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation; either version 2 of the License, or
## (at your option) any later version.
##
## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with this program; if not, write to the Free Software
## Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

from numpy import *
from scipy.ndimage import convolve
__all__ = ['bweuler']

def bweuler(bw,n=8):
    """
    Calculates the Euler number of a binary image

    eul=bweuler(bw,n) calculates the Euler number @var{eul} of a binary
    image bw, which is a scalar whose value is the total number of
    objects in an image minus the number of holes.

    n can have one of two values:
      4: bweuler will use 4-connected neighbourhood definition.
      8: bweuler will use 8-connected neighbourhood definition. This is the default value.

    This function uses Bit Quads as described in "Digital Image
    Processing" to calculate euler number.

    References:
    W. K. Pratt, "Digital Image Processing", 3rd Edition, pp 593-595
    """
    if n==8:
        lut=array([0,.25,.25,0,.25,0,-.5,-.25,.25,-.5,0,-.25,0,-.25,-.25,0])
    elif n==4:
        lut=array([0,.25,.25,0,.25,0,.5,-.25,.25,.5,0,-.25,0,-.25,-.25,0])
    else:
        raise Exception("bweuler: n can only be 4 or 8.");
      
    A=applylut(bw,lut)
    return A.sum()


def applylut(bw,lut):
    """
    A= applylut (bw,lut)
    Uses lookup tables to perform a neighbour operation on binary images.

    A = applylut(BW,LUT) returns the result of a neighbour operation
    using the lookup table LUT which can be created by makelut.

    It first computes a matrix with the index of each element in the
    lookup table. To do this, it convolves the original matrix with a
    matrix which assigns each of the neighbours a bit in the resulting
    index. Then LUT is accessed to compute the result.
    """
    nq=log2(lut.size)
    n=sqrt(nq)
    if floor(n)!=n:
        raise Exception("applylut: LUT length is not as expected. Use makelut to create it.")
    w=reshape(2**(nq-1-arange(nq)),(n,n))
    A=lut[convolve(asarray(bw,uint8),w)]
    return A

# vim: set ts=4 sts=4 sw=4 expandtab smartindent:
