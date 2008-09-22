from string import upper
from ..image import Image
from ..preprocess import preprocessimage
from numpy import *

from edgefeatures import edgefeatures
from texture import haralickfeatures
from imgskelfeats import imgskelfeatures
from noffeatures import noffeatures
from imgfeatures import imgfeatures, imgfeaturesdna
from hullfeatures import hullfeatures, hullsizefeatures
from zernike import zernike

__all__ = ['computefeatures','featurenames']

def _featsfor(featset):
    ufeatset=upper(featset)
    if ufeatset == 'SLF13' or ufeatset == 'SLF7DNA':
        return ['skl','nof','img','hul','zer','har','edg']
    if ufeatset == 'MCELL':
        return ['har','img','edg','skl']
    return [featset]

def computefeatures(img,featsets):
    '''
    Compute features on img.

    featsets can be either a list of feature groups (currently recognised:
        'skl' 'nof' 'img' 'hul' 'zer' 'har' 'edg') or a feature set name
        (currently recognised 'SLF7dna', 'mcell')

    img can be a list of images. In this case, a list of feature vectors will be returned.
        Also, in this case, imgs will be unload after feature calculation.
    '''
    if type(featsets) == str:
        featsets = _featsfor(featsets)
    if type(img) == list:
        features=[]
        for i in img:
            f=computefeatures(i,featsets)
            features.append(f)
            i.unload()
        return features
    scale=img.scale
    if scale is None:
        scale = .23
    preprocessimage(img,1,{})
    features=array([])
    procprotein=img.channeldata[Image.procprotein_channel]
    resprotein=img.channeldata[Image.residualprotein_channel]
    procdna=img.channeldata.get(Image.procdna_channel)
    for F in featsets:
        if F == 'edg':
            feats=edgefeatures(procprotein)
        elif F == 'har':
            feats=haralickfeatures(procprotein)
            feats=feats.mean(0)
        elif F == 'hul':
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
        elif F == 'skl':
            feats=imgskelfeatures(procprotein)
        elif F == 'zer':
            feats=zernike(procprotein,12,34.5,scale)
        else:
            raise Exception('Unknown feature set: %s' % F)
        features = r_[features,feats]
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
        else:
            raise Exception('Unknown feature set: %s' % F)
    return names

# vim: set ts=4 sts=4 sw=4 expandtab smartindent:
