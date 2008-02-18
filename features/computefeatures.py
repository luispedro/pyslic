from sys import path
path.append('../')
from string import upper
from image import Image
from preprocess import preprocess
from numpy import *

from edgefeatures import edgefeatures
from texture import haralickfeatures
from imgskelfeats import imgskelfeatures
from noffeatures import noffeatures
from imgfeatures import imgfeatures
from hullfeatures import hullfeatures
from zernike import zernike

__all__ = ['computefeatures']

def _featsfor(featset):
    featset=upper(featset)
    if featset == 'SLF13' or featset == 'SLF7DNA':
        return ['skl','nof','img','hul','zer','har','edg']

def computefeatures(img,featsets):
    if type(featsets) == str:
        featsets = _featsfor(featsets)
    if type(img) == list:
        features=[]
        for i in img:
            f=computefeatures(i,featsets)
            features.append(f)
        return features
    preprocess(img,1,{})
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
        elif F == 'img':
            feats=imgfeatures(procprotein,procdna)
        elif F == 'mor':
            feats=morphologicalfeatures(procprotein)
        elif F == 'nof':
            feats=noffeatures(procprotein,resprotein)
        elif F == 'skl':
            feats=imgskelfeatures(procprotein)
        elif F == 'zer':
            feats=zernike(procprotein,12,34.5)
        else:
            raise Exception('Unkown feature set: %s' % F)
        features = r_[features,feats]
    return features

# vim: set ts=4 sts=4 sw=4 expandtab smartindent:
