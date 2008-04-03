from string import upper
from ..image import Image
from ..preprocess import preprocessimage
from numpy import *

from edgefeatures import edgefeatures
from texture import haralickfeatures
from imgskelfeats import imgskelfeatures
from noffeatures import noffeatures
from imgfeatures import imgfeatures
from hullfeatures import hullfeatures, hullsizefeatures
from zernike import zernike

__all__ = ['computefeatures']

def _featsfor(featset):
    ufeatset=upper(featset)
    if ufeatset == 'SLF13' or ufeatset == 'SLF7DNA':
        return ['skl','nof','img','hul','zer','har','edg']
    return [featset]

def computefeatures(img,featsets):
    '''
    Compute features on img.

    featsets can be either a list of feature groups (currently recognised:
        'skl' 'nof' 'img' 'hul' 'zer' 'har' 'edg') or a feature set name
        (currently recognised 'SLF7dna')

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
            feats=imgfeatures(procprotein,procdna)
        elif F == 'mor':
            feats=morphologicalfeatures(procprotein)
        elif F == 'nof':
            feats=noffeatures(procprotein,resprotein)
        elif F == 'skl':
            feats=imgskelfeatures(procprotein)
        elif F == 'zer':
            feats=zernike(procprotein,12,34.5,scale)
        else:
            raise Exception('Unkown feature set: %s' % F)
        features = r_[features,feats]
    return features

# vim: set ts=4 sts=4 sw=4 expandtab smartindent:
