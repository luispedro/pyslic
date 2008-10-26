# -*- coding: utf-8 -*-
# Copyright (C) 2008  Murphy Lab
# Carnegie Mellon University
# 
# Written by Luís Pedro Coelho <lpc@cmu.edu>
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
import numpy
from scipy import ndimage

__all__ = ['pftas','tas']
_M = numpy.array([
    [1, 1,1],
    [1,10,1],
    [1, 1,1]
    ])
_bins = numpy.arange(11)

def _tas(img,thresh,margin):
    def _ctas(img):
        V = ndimage.convolve(img,_M)
        values,_ = numpy.histogram(V,bins=_bins)
        values = values[:9]
        return values/values.sum()

    mu = ((img > thresh)*img).sum() / (img > thresh).sum()
    bimg = (img > mu - margin) * (img < mu + margin)
    return numpy.concatenate((_ctas(bimg),_ctas(~bimg)))

def tas(img):
    '''
    values = tas(img)

     This algorithm was presented by Hamilton et al.
    in "Fast automated cell phenotype image classification"
    (http://www.biomedcentral.com/1471-2105/8/110)

    See also pftas() for a variation without any hardcoded parameters.
    '''
    return _tas(img,30,30)
tas.names = [( 'tas_%s' % i) for i in xrange(9)] + \
            [('ntas_%s' % i) for i in xrange(9)]

def pftas(img):
    '''
    values = pftas(img)

     This algorithm was presented by Hamilton et al.
    in "Fast automated cell phenotype image classification"
    (http://www.biomedcentral.com/1471-2105/8/110)

     The current implementation is an adapted version which
    is free of parameters. The thresholding is done beforehand
    (automatically), the margin around the mean of pixels to be
    included is the standard deviation of the pixel values
    and not fixed to 30, as before.

     Use tas() to get the original version of the features.
    Also do NOT run the preprocessing code on images on which the original
    are to be calculated on.
    '''
    T=0
    std=img[img>T].std()
    return _tas(img,T,std)

pftas.names = [( 'pftas_%s' % i) for i in xrange(9)] + \
              [('npftas_%s' % i) for i in xrange(9)]

# vim: set ts=4 sts=4 sw=4 expandtab smartindent:
