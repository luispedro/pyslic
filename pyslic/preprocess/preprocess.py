# -*- coding: utf-8 -*-
# vim: set ts=4 sts=4 sw=4 expandtab smartindent:
# Copyright (C) 2008-2010 Murphy Lab
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
import numpy as np
from ..image import Image
from ..imageprocessing.thresholding import rc
from mahotas.morph import majority_filter
from mahotas.histogram import fullhistogram
from mahotas.bbox import bbox
from mahotas.stretch import stretch
from scipy import ndimage
from warnings import warn
fn = np

__all__ = ['preprocessimg', 'precomputestats', 'bgsub']

def precomputestats(image):
    image.lazy_load()
    image.temp['bgsubprotein'] = bgsub(image.channeldata['protein'].copy())
    if 'dna' in image.channeldata:
        image.temp['bgsubdna'] = bgsub(image.channeldata['dna'].copy())
    if image.regions is not None:
        image.temp['region_ids'] = ndimage.find_objects(image.regions)


def preprocessimage(image, regionid=None, crop=True, options = {}):
    """
    Preprocess the image

    image should be an Image object
    regionid should be an integer, which indexes into image.region

    options is a dictionary. The following options are accepted:

        + 'bgsub.way': One of
            - 'ml' (default) 
            - 'mb' and controls whether the region is selected prior to preprocessing
         + '3d.mode': One of
            - 'perslice' (default): process each slice separately
         + 'threshold.algorithm': One of
            - 'rc': Riddlar-Calvard
            - 'mean':  Image mean
    """
    def preprocessimg(img):
        if len(img.shape) > 2:
            assert len(img.shape) == 3, "Cannot handle images of more than 3 dimensions."
            if options.get('3d.mode','perslice') == 'perslice':
                nr_slices=img.shape[0]
                out_proc=img.copy()
                out_res=img.copy()
                for z in xrange(nr_slices):
                    proc,res=preprocessimg(img[z])
                    out_proc[z]=proc
                    out_res[z]=res
                return out_proc,out_res
            else:
                raise Exception('pyslic.preprocessimg: Do not know how to handle 3d.mode: %s' % options['3d.mode'])
        if do_bgsub:
            regions = image.regions
            img = img.copy()
            if regions is not None:
                if options.get('bgsub.way','ml') == 'ml':
                    img *= (regions == regionid)
                    img = bgsub(img, options)
                else:
                    img = bgsub(img, options)
                    img *= (cropimg == regionid)
            else:
                if regionid:
                    warn('Selecting a region different from 1 for an image without region information')
                img = bgsub(img, options)
        imgscaled = stretch(img, 255)
        T = thresholdfor(imgscaled,options)
        mask = (imgscaled > T)
        mask = majority_filter(mask)
        residual = img.copy()
        img *= mask
        residual *= ~mask
        return img,residual
    image.lazy_load()

    protein = image.channeldata['protein']
    dna = image.channeldata.get('dna')
    do_bgsub = True
    if 'bgsubprotein' in image.temp:
        protein = image.temp['bgsubprotein']
        do_bgsub = False
    if regionid is not None and 'region_ids' in image.temp:
        location = image.temp['region_ids'][regionid - 1]
        if location is None:
            image.channeldata['procprotein'] = \
                image.channeldata['resprotein'] = \
                image.channeldata['procdna'] = np.zeros((0,0), dtype=protein.dtype)
            return
        protein = protein[location]
        if dna is not None:
            dna = dna[location]
    image.channeldata['procprotein'],image.channeldata['resprotein'] = preprocessimg(protein)
    if dna is not None:
        image.channeldata['procdna'],_ = preprocessimg(dna)

    if crop:
        fullimage = (image.channeldata['procprotein'] > 0) | (image.channeldata['resprotein'] >0)
        if 'dna' in image.channeldata:
            fullimage |= (image.channeldata['procdna'] > 0)

        min1,max1,min2,max2 = bbox(fullimage)
        border = 2
        min1 = max(0, min1 - border)
        min2 = max(0, min2 - border)
        max1 += border
        max2 += border

        image.channeldata['procprotein'] = image.channeldata['procprotein'][min1:max1,min2:max2]
        image.channeldata['resprotein'] = image.channeldata['resprotein'][min1:max1,min2:max2]
        if 'dna' in image.channeldata:
            image.channeldata['procdna'] = image.channeldata['procdna'][min1:max1,min2:max2]


def thresholdfor(img,options = {}):
    type = options.get('threshold.algorithm','rc')
    if type == 'rc':
        return rc(img,ignore_zeros=True)
    elif type == 'mean':
        return img.mean()
    else:
        raise KeyError('Threshold option not recognised (%s).' % type)

def bgsub(img,options = {}):
    '''
    bgsub(img,options = None)

    Background subtract img (which must be a numpy-type array).

    Changes are done inplace and the img is returned. Use the following idiom for operating on a copy:

    B = bgsub(A.copy(),options)
    '''
    type = options.get('bgsub.type','lowcommon')
    if type == 'nobgsub':
        return img
    elif type == 'lowcommon':
        M = np.round(img.mean())-1
        if M <= 0:
            T = 0
        else:
            hist = fullhistogram(img)
            T = np.argmax(hist[:M])
        if T > 0:
            img -= np.minimum(img,T)
        return img
    else:
        raise KeyError('Background subtraction option not recognised (%s).' % type)

