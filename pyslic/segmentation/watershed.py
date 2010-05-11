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

from __future__ import division, with_statement
import numpy as np
from ..image import loadedimage
from ..imageprocessing.thresholding import threshold
import mahotas
import pymorph
from scipy import ndimage

__all__ = ['watershed']

def watershed_segment(img, mode='direct', thresholding=None, min_obj_size=None, **kwargs):
    '''
    segment_watershed(img, mode='direct', thresholding=None, min_obj_size=None, **kwargs)

    Segment using traditional watershed

    Parameters
    ----------

        * img: a pyslic.Image. The algorithm operates on the dna channel.
        * mode: 'direct' or 'gradient': whether to use the image or the gradient of the image
        * thresholding: how to threshold the smoothed image (default: None, no thresholding)
        * min_obj_size: minimum object size. This is slightly different than post-filtering for minimum 
            object size as it fill those holes with a second watershed pass as opposed to having
            an image with holes
        * smoothing: whether to smooth (default: True)
        * smooth_gamma: Size of Gaussian for blurring, in pixels (default: 12)
    '''
    assert mode in ('direct','gradient'), "segment_watershed: mode '%s' not understood" % mode
    with loadedimage(img):
        dna = img.get('dna')
        if kwargs.get('smoothing',True):
            dnaf = ndimage.gaussian_filter(dna, kwargs.get('smooth_gamma',12))
        else:
            dnaf = dna
        rmax = pymorph.regmax(dnaf)
        rmax_L,_ = ndimage.label(rmax)
        if mode == 'direct':
            watershed_img = dna.max()-dna
        elif mode == 'gradient':
            dnag = pymorph.gradm(dna)
            watershed_img = dnag.max()-dnag
        water = mahotas.cwatershed(watershed_img,rmax_L)
        if thresholding is not None:
            T = threshold(dnaf,thresholding)
            water *= (dnaf >= T)
        if min_obj_size is not None:
            oid = 1
            while oid <= water.max():
                if (water == oid).sum() < min_obj_size:
                    water[water == oid] =0
                    water[water > oid] -= 1
                else:
                    oid += 1
        return water

# vim: set ts=4 sts=4 sw=4 expandtab smartindent:
