# Copyright (C) 2006  Murphy Lab
# Carnegie Mellon University
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

# 19 Dec 98 - M.V. Boland

from __future__ import division
from math import *
from numpy import *
from scipy.ndimage import *

__all__ = ['zernike']

def _factorial(N):
    res = 1.
    for i in xrange(1,N+1):
        res *= i
    return res

def _polar(r,theta):
    x = r * cos(theta)
    y = r * sin(theta)
    return 1*x+1j*y

def Znl(n,l,X,Y,P):
    v = 0.+0.j
    for x,y,p in zip(X,Y,P):
        Vnl = 0.
        for m in xrange( (n-l)//2 + 1 ):
              Vnl += (-1.)**m * _factorial(n-m) /  \
            ( _factorial(m) * _factorial((n - 2*m + l) // 2) * _factorial((n - 2*m - l) // 2) ) * \
            ( sqrt(x*x + y*y)**(n - 2*m) * _polar(1.0, l*atan2(y,x)) )
        v += p * conjugate(Vnl)

    v *= (n+1)/pi
    return v 


def zernike(img,D,radius):
    """
     ZVALUES = zernike(img,D,radius) Zernike moments through degree D 
     ML_ZERNIKE(I,D,R),
     Returns a vector of Zernike moments through degree D for the
     image I, and the names of those moments in cell array znames. 
     R is used as the maximum radius for the Zernike polynomials.

     For use as features, it is desirable to take the 
     magnitude of the Zernike moments (i.e. abs(zvalues))

     Reference: Teague, MR. (1980). Image Analysis via the General
       Theory of Moments.  J. Opt. Soc. Am. 70(8):920-930.
    """
    znames = []
    zvalues = []

# Find all non-zero pixel coordinates and values
    X,Y = where(img > 0)
    P=img[X,Y].ravel()

# Normalize the coordinates to the center of mass and normalize
#  pixel distances using the maximum radius argument (radius)
    cofx,cofy = center_of_mass(img)
    Xn = double(X-cofx)/radius
    Yn = double(Y-cofy)/radius
    Xn=Xn.ravel()
    Yn=Yn.ravel()


# Find all pixels of distance <= 1.0 to center
    k = (sqrt(Xn**2 + Yn**2) <= 1.)
    frac_center = array(P[k],double)/img.sum()

    for n in xrange(D+1):
        for l in xrange(n+1):
            if (n-l)%2 == 0:
                znames.append('Z_#i,%s%s' % (n, l))
                z= Znl(n,l, Xn[k], Yn[k], frac_center.ravel())
                zvalues.append(abs(z))
    return zvalues

# vim: set ts=4 sts=4 sw=4 expandtab smartindent:
