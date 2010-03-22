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
from numpy import *
from scipy.ndimage import median_filter, label, center_of_mass, gaussian_filter
from ..imageprocessing.thresholding import otsu, murphy_rc

__all__ = ['labelnuclei','nucleicof']

def labelnuclei(dnaimg,**options):
    '''
    labeled,N = labeled_nuclei(dnaimg, **options)

    N equals the number of nuclei
    labeled is of the same shape as dnaimg and contains the labels of the nuclei
    '''
    sigma=options.get('sigma',10)
    dnaimg=gaussian_filter(dnaimg,sigma)
    thresholding=options.get('thresholding','otsu')
    if thresholding == 'otsu':
        T=otsu(dnaimg)
    elif thresholding == 'murphy_rc':
        T=murphy_rc(dnaimg)
    else:
        raise AttributeError, "Unknown thresholding method (%s)" % thresholding
    return label(dnaimg > T)

def nucleicof(dnaimg,options=None):
    '''
    Returns a set of nuclear centres.
    '''
    labeled,N=labelnuclei(dnaimg)
    cofs=center_of_mass(dnaimg,labeled,range(1,N+1))
    return cofs

# vim: set ts=4 sts=4 sw=4 expandtab smartindent:
