# -*- coding: utf-8 -*-
# Copyright (C) 2008  Murphy Lab
# Carnegie Mellon University
# 
# Original (Matlab) version by M.V. Boland (10 Aug 98)
# Python version by Luís Pedro Coelho <lpc@cmu.edu>
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

from numpy import *
from scipy.ndimage import center_of_mass

__all__ = ['imgmoments','imgcentmoments']
def imgmonents(img,x,y):
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
    img = double(img)
    N,M=img.shape
    img *= arange(M)**y
    img=img.T
    img *= arange(N)**x
    return img.sum()

def imgcentmoments(img,x,y,cofy=None,cofx=None):
    """
    M_xy = imgcentmoments(img,x,y, cofy=None, cofx=None)

    @param cofy and cofx are optional and computed from the image if not given
    """
    if cofy is None or cofx is None:
        print 'calling center_of_mass'
        cofy, cofx = center_of_mass(img)
    try:
        from scipy import weave
        from scipy.weave import converters
        max1,max2=img.shape
        img=asarray(img,uint8)
        vals=zeros(1)
        code = '''
#line 67 "imgmoments.py"
        double res = 0.;
        for (int y_index = 0; y_index != max1; ++y_index) {
            for (int x_index = 0; x_index != max2; ++x_index) {
                res += img(y_index,x_index) * std::pow(cofy - y_index,y) * std::pow(cofx - x_index,x);
            }
        }
        vals(0) = res;
        '''
        weave.inline(
                code,
                ['max1','max2','img','cofy','cofx','y','x','vals'],
                type_converters=converters.blitz)
        return vals[0]
    except:
        import warnings
        warnings.warn('scipy.weave failed. Resorting to (slow) Python code')
        img = double(img)
        N,M=img.shape
        img *= (arange(M)-cofx)**x
        img = img.T
        img *= (arange(N)-cofy)**y
        return img.sum()

# vim: set ts=4 sts=4 sw=4 expandtab smartindent:
