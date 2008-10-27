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
from string import upper
from ..image import Image
from ..preprocess import preprocessimage

from edgefeatures import edgefeatures
from texture import haralickfeatures
from imgskelfeats import imgskelfeatures
from noffeatures import noffeatures
from imgfeatures import imgfeatures, imgfeaturesdna
from hullfeatures import hullfeatures, hullsizefeatures
from zernike import zernike
from tas import tas, pftas

__all__ = ['computefeatures','featurenames']

def _featsfor(featset):
    ufeatset=upper(featset)
    if ufeatset == 'ALL':
        return ['skl','nof','img','hul','zer','har','edg','pftas']
    if ufeatset == 'SLF13' or ufeatset == 'SLF7DNA':
        return ['skl','nof','img','hul','zer','har','edg']
    if ufeatset == 'MCELL':
        return ['har','img','edg','skl']
    return [featset]

_Default_Haralick_Scale = 1.15
_Default_Haralick_Bins = 32

def computefeatures(img,featsets,progress=None,**kwargs):
    '''

    features = computefeatures(img,featsets,progress=None)

    Compute features defined by featsets on image img.

    featsets can be either a list of feature groups.
    Feature groups:
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

    img can be a list of images. In this case, a two-dimensional feature vector will be returned, where
    f[i,j] is the j-th feature of the i-th image. Also, in this case, imgs will be unload after feature calculation.

    @param progress: if progress is not None, then it should be an integer. Every progress images, an output message will be printed.
    '''
    if type(featsets) == str:
        featsets = _featsfor(featsets)
    if type(img) == list:
        features=[]
        for i,im in enumerate(img):
            f=computefeatures(im,featsets)
            features.append(f)
            im.unload()
            if progress is not None and (i % progress) == 0:
                print 'Processed %s images...' % i
        return numpy.array(features)
    scale=img.scale
    if scale is None:
        scale = .23
    preprocessimage(img,1,{})
    features=numpy.array([])
    protein=img.channeldata[Image.protein_channel]
    procprotein=img.channeldata[Image.procprotein_channel]
    resprotein=img.channeldata[Image.residualprotein_channel]
    procdna=img.channeldata.get(Image.procdna_channel)
    for F in featsets:
        if F in ['edg','edge']:
            feats=edgefeatures(procprotein)
        elif F == 'har':
            img=procprotein
            har_scale = kwargs.get('haralick.scale',_Default_Haralick_Scale)
            if scale != har_scale:
                img = img.copy()
                img = ndimage.zoom(img, scale/_Default_Haralick_Scale)
            bins = kwargs.get('haralick.bins',_Default_Haralick_Bins)
            if bins != 256:
                img = numpy.array(img.astype(float) * bins / (img.max()-img.min()),numpy.uint8)
            feats=haralickfeatures(procprotein)
            feats=feats.mean(0)
        elif F == 'har3d':
            feats=haralickfeatures(procprotein)
            feats=numpy.r_[feats.mean(0),feats.ptp(0)]
        elif F in ['hul', 'hull']:
            feats=hullfeatures(procprotein)
        elif F == 'hullsize':
            feats=hullsizefeatures(procprotein)
        elif F == 'hullsizedna':
            feats=hullsizefeatures(procdna)
        elif F == 'img':
            feats=imgfeaturesdna(procprotein,procdna)
        elif F == 'mor':
            feats=morphologicalfeatures(procprotein)
        elif F == 'nof':
            feats=noffeatures(procprotein,resprotein)
        elif F in ['skl', 'skel']:
            feats=imgskelfeatures(procprotein)
        elif F == 'zer':
            feats=zernike(procprotein,12,34.5,scale)
        elif F == 'tas':
            feats=tas(protein)
        elif F == 'pftas':
            feats=pftas(procprotein)
        else:
            raise Exception('Unknown feature set: %s' % F)
        features = numpy.r_[features,feats]
    return features

def featurenames(featsets):
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
            feats=(procprotein)
            names.extend(zernike.names)
        elif F == 'tas':
            names.extend(tas.names)
        elif F == 'pftas':
            names.extend(pftas.names)
        else:
            raise Exception('Unknown feature set: %s' % F)
    return names

# vim: set ts=4 sts=4 sw=4 expandtab smartindent:
