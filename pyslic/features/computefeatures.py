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
import numpy as np
import numpy
from scipy import ndimage
from string import upper
import re
from ..image import Image
from ..preprocess import preprocessimage, precomputestats

from edgefeatures import edgefeatures
from texture import haralickfeatures
from imgskelfeats import imgskelfeatures
from noffeatures import noffeatures
from imgfeatures import imgfeatures, imgfeaturesdna
from hullfeatures import hullfeatures, hullsizefeatures
from zernike import zernike, znames
from tas import tas, pftas
from .overlap import overlapfeatures
from mahotas.lbp import lbp

__all__ = ['computefeatures','featurenames']

def _featsfor(featset):
    ufeatset=upper(featset)
    if ufeatset == 'ALL':
        return ['skl','nof','img','hul','zer','har','edg','pftas']
    if ufeatset == 'SLF13' or ufeatset == 'SLF7DNA':
        return ['skl','nof','img','hul','zer','har','edg']
    if ufeatset == 'MCELL':
        return ['har','obj-field','edg','skl']
    if ufeatset == 'SLF33':
        return ['har','har1','har2','har3','har4','har5','har6','obj-field','edg','skl','nof','pftas']
    if ufeatset == 'SLF34':
        return ['har','har1','har2','har3','har4','har5','har6','obj-field-dna','edg','skl','nof','pftas','overlap']
    if ufeatset == 'FIELD+':
        return ['har','har1','har2','har3','har4','har5','har6','obj-field','edg','skl','nof','pftas']
    if ufeatset == 'FIELD-DNA+':
        return ['har','har1','har2','har3','har4','har5','har6','obj-field-dna','edg','skl','nof','pftas','overlap']
    return [featset]

_Default_Scale = .23
_Default_Haralick_Scale = 1.15
_Default_Haralick_Bins = 32
_Min_image_size = 100

def computefeatures(img, featsets, progress=None, preprocessing=True, **kwargs):
    '''
    features = computefeatures(img,featsets,progress=None,**kwargs)

    Compute features defined by featsets on image img.


    featsets can be either a list of feature groups.

    Feature groups:
    --------------
        + 'har': 13 Haralick features [in matslic]
        + 'har3d': 26 Haralick features (actually this works in 2D as well)
        + 'skl': Skeleton features [in matslic] (syn: 'skel')
        + 'img': Image features [in matslic]
        + 'hul': Hull features [in matslic] (syn: 'hull')
        + 'edg': Edge features [in matslic] (syn: 'edge')
        + 'zer': Zernike moments [in matslic]
        + 'tas' : Threshold Adjacency Statistics
        + 'pftas' : Parameter-free Threshold Adjacency Statistics
    Feature set names:
        + 'SLF7dna'
        + 'mcell': field level features
        + 'field+': all field level features. The exact exact number of
                    features computed is liable to change (increase) in
                    newer versions

    img can be a list of images. In this case, a two-dimensional feature vector will be returned, where
    f[i,j] is the j-th feature of the i-th image. Also, in this case, imgs will be unload after feature calculation.

    Parameters
    ----------
        * *progress*: if progress is not None, then it should be an integer.
                Every *progress* images, an output message will be printed.
        * *preprocessing*: Whether to do preprocessing (default: True).
                If False and img.channeldata[procprotein|procdna] are empty, they are filled in
                using img.channeldata[protein|dna] respectively.
        * *options*: currently passed through to pyslic.preprocessimage
    '''
    if type(featsets) == str:
        featsets = _featsfor(featsets)
    if type(img) == list:
        features=[]
        for i,im in enumerate(img):
            f = computefeatures(im,featsets,progress=None,**kwargs)
            features.append(f)
            im.unload()
            if progress is not None and (i % progress) == 0:
                print 'Processed %s images...' % i
        return numpy.array(features)
    regions = img.regions
    if regions is not None and regions.max() > 1 and 'region' not in kwargs:
        precomputestats(img)
        return numpy.array([
                    computefeatures(img, featsets, progress=progress, region=r, **kwargs)
                    for r in xrange(1,regions.max()+1)])
    scale = img.scale
    if scale is None:
        scale = _Default_Scale
    if preprocessing:
        preprocessimage(img, kwargs.get('region',1), options=kwargs.get('options',{}))
    else:
        if not img.channeldata['procprotein']:
            img.channeldata['procprotein'] = img.get('protein')
            img.channeldata['resprotein'] = 0*img.get('protein')
        if not img.channeldata['procdna']:
            img.channeldata['procdna'] = img.get('dna')
    features = numpy.array([])
    protein = img.get('protein')
    procprotein = img.get('procprotein')
    resprotein = img.get('resprotein')
    dna = img.channeldata.get('dna')
    procdna = img.channeldata.get('procdna')
    if procprotein.size < _Min_image_size:
        return np.array([np.nan for i in xrange(90)])
    lbppat = re.compile(r'lbp\(([0-9]+), ?([0-9]+)\)')
    for F in featsets:
        if F in ['edg','edge']:
            feats = edgefeatures(procprotein)
        elif F == 'raw-har':
            feats = haralickfeatures(procprotein).mean(0)
        elif F[:3] == 'har':
            img = procprotein
            har_scale = kwargs.get('haralick.scale',_Default_Haralick_Scale)
            if F == 'har' and scale != har_scale:
                img = img.copy()
                img = ndimage.zoom(img, scale/_Default_Haralick_Scale)
            if len(F) > 3:
                rate = int(F[3])
                if rate != 1:
                    C = np.ones((rate,rate))
                    img = np.array(img,np.uint16)
                    img = ndimage.convolve(img,C)
                    img = img[::rate,::rate]
            if not img.size:
                feats = np.zeros(13)
            else:
                bins = kwargs.get('haralick.bins',_Default_Haralick_Bins)
                if bins != 256:
                    min = img.min()
                    max = img.max()
                    ptp = max - min
                    if ptp:
                        img = np.array((img-min).astype(float) * bins/ptp, np.uint8)
                feats = haralickfeatures(img)
                feats = feats.mean(0)
        elif F == 'har3d':
            feats = haralickfeatures(procprotein)
            feats = numpy.r_[feats.mean(0),feats.ptp(0)]
        elif F in ['hul', 'hull']:
            feats = hullfeatures(procprotein)
        elif F == 'hullsize':
            feats = hullsizefeatures(procprotein)
        elif F == 'hullsizedna':
            feats = hullsizefeatures(procdna)
        elif F == 'img':
            feats = imgfeaturesdna(procprotein, procdna)
        elif F in ('obj-field', 'obj-field-dna'):
            feats = imgfeaturesdna(procprotein, procdna, isfield=True)
        elif F == 'mor':
            feats = morphologicalfeatures(procprotein)
        elif F == 'nof':
            feats = noffeatures(procprotein,resprotein)
        elif F in ['skl', 'skel']:
            feats = imgskelfeatures(procprotein)
        elif F == 'zer':
            feats = zernike(procprotein,12,34.5,scale)
        elif F == 'tas':
            feats = tas(protein)
        elif F == 'pftas':
            feats = pftas(procprotein)
        elif F == 'overlap':
            feats = overlapfeatures(protein, dna, procprotein, procdna)
        elif lbppat.match(F):
            radius,points = lbppat.match(F).groups()
            feats = lbp(protein, int(radius), int(points))
        else:
            raise Exception('Unknown feature set: %s' % F)
        features = numpy.r_[features,feats]
    return features

def featurenames(featsets):
    '''
    names = featurenames(featsets)

    Returns a list of feature names. The argument has the same
    meaning as the argument to computefeatures.
    '''
    featsets = _featsfor(featsets)
    names=[]
    for F in featsets:
        if F == 'edg':
            names.extend(edgefeatures.names)
        elif F == 'har':
            names.extend(haralickfeatures.names)
        elif F == 'hul':
            names.extend(hullfeatures.names)
        elif F == 'hullsize':
            names.extend(hullsizefeatures.names)
        elif F == 'hullsizedna':
            names.extend(hullsizefeatures.names)
        elif F == 'img':
            names.extend(imgfeaturesdna.names)
        elif F == 'imgnodna':
            names.extend(imgfeatures.names)
        elif F == 'imgdna':
            names.extend(imgfeaturesdna.names)
        elif F == 'mor':
            names.extend(morphologicalfeatures.names)
        elif F == 'nof':
            names.extend(noffeatures.names)
        elif F == 'skl':
            names.extend(imgskelfeatures.names)
        elif F == 'zer':
            names.extend(znames(12,34.5))
        elif F == 'tas':
            names.extend(tas.names)
        elif F == 'pftas':
            names.extend(pftas.names)
        else:
            raise Exception('Unknown feature set: %s' % F)
    return names

# vim: set ts=4 sts=4 sw=4 expandtab smartindent:
