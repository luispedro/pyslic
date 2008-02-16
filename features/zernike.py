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
from scipy import weave
from scipy.weave import converters

__all__ = ['zernike']

def _polar(r,theta):
    x = r * cos(theta)
    y = r * sin(theta)
    return 1*x+1j*y

def Znl(n,l,X,Y,P):
    v = 0.+0.j
    factorialtable=array([1,1,2,6,24,120,720,5040,40320,362880,3628800,39916800,479001600])
    try:
        Nelems=len(X)
        v=array([v]) # This is necessary for the C++ code to see and update v correctly
        code='''
#line 45 "zernike.py"
        using std::pow;
        using std::atan2;
        using std::polar;
        using std::conj;
        using std::complex;
        complex<double> Vnl = 0.0;
        v(0)=0;
        for (int i = 0; i != Nelems; ++i) {
            double x=X(i);
            double y=Y(i);
            double p=P(i);
            Vnl = 0.;
            for(int m = 0; m <= (n-l)/2; m++) {
                double f = (m & 1) ? -1 : 1;
                Vnl += f * factorialtable(int(n-m)) /
                       ( factorialtable(m) * factorialtable((n - 2*m + l) / 2) * factorialtable((n - 2*m - l) / 2) ) *
                       ( pow( sqrt(x*x + y*y), (double)(n - 2*m)) ) *
                       polar(1.0, l*atan2(y,x)) ;
            }
            v(0) += p * conj(Vnl);
        }
        '''
        weave.inline(code,
            ['factorialtable','X','Y','P','v','n','l','Nelems'],
            type_converters=converters.blitz,
            compiler = 'gcc',
            headers=['<complex>'])
        v=v[0]
    except:
        for x,y,p in zip(X,Y,P):
            Vnl = 0.
            for m in xrange( (n-l)//2 + 1 ):
                  Vnl += (-1.)**m * factorialtable[n-m] /  \
                ( factorialtable[m] * factorialtable[(n - 2*m + l) // 2] * factorialtable[(n - 2*m - l) // 2] ) * \
                ( sqrt(x*x + y*y)**(n - 2*m) * _polar(1.0, l*atan2(y,x)) )
            v += p * conjugate(Vnl)
    v *= (n+1)/pi
    return v 


def zernike(img,D,radius,scale=.23):
    """
     zvalues = zernike(img,D,radius) zernike moments through degree D

     Returns a vector of absolute Zernike moments through degree D for the
     image I, and the names of those moments in cell array znames. 
     radius is used as the maximum radius for the Zernike polynomials.

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
    Xn = double(X-cofx)/radius*scale
    Yn = double(Y-cofy)/radius*scale
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
