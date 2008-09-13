from __future__ import division
from ..image import Image
from numpy import *
from ..imageprocessing.bweuler import bweuler
from ..imageprocessing.bbox import bbox
from scipy import ndimage
from scipy.ndimage import *
from .hullfeatures import hullfeatures
from .mmthin import mmthin
from ..imageprocessing.convexhull import convexhull
from .imgskelfeats import find_branch_points
import warnings

def objectfeatures(img,region=None):
    '''
    values=objectfeatures(img,region=None)

    This implements the object features described in
    "Object Type Recognition for Automated Analysis of Protein Subcellular Location"
    by Ting Zhao, Meel Velliste, Michael V. Boland, and Robert F. Murphy
    in IEEE Transaction on Image Processing
    '''

    protimg = img.channeldata[Image.procprotein_channel]
    dnaimg = img.channeldata.get(Image.procdna_channel,None)
    if region is not None:
        if img.regions is None:
            assert region == 1, 'objectfeatures called with region != 1, but img.regions is None!'
        else:
            protimg = protimg * (img.regions == region)
            dnaimg = dnaimg * (img.regions == region)

    labeled,N = ndimage.label(protimg,ones((3,3)))
    objects = xrange(1,N+1)

    sofs = zeros((len(objects),11))
    if dnaimg is not None:
        dnacofy,dnacofx=center_of_mass(dnaimg)
        bindna=(dnaimg > 0)
    for obji,obj in enumerate(objects):
        binobj = (labeled == obj)
        min1,max1,min2,max2=bbox(binobj)
        if min1 > 0: min1 -= 1
        if min2 > 0: min2 -= 1
        binobjc = binobj[min1:max1+1,min2:max2+1] # leave a small margin for bweuler()
        protobj = protimg[min1:max1+1,min2:max2+1]
        objimg = protimg * binobj
        cofy,cofx = center_of_mass(objimg)
        objskel = mmthin(binobjc)
        binskel = (objskel > 0)
        objhull=convexhull(binobjc)
        no_of_branch_points = find_branch_points(objskel).sum()
        hfeats=hullfeatures(binobjc,objhull)

        sofs[obji,0] = binobjc.sum()
        if dnaimg is not None:
            sofs[obji,1] = sqrt((cofy-dnacofy)**2+(cofx-dnacofx)**2)
            sofs[obji,2] = (binobj*bindna).sum()/sofs[obji,0]
        sofs[obji,3] = hfeats[2]
        sofs[obji,4] = bweuler(binobjc)
        sofs[obji,5] = hfeats[1]
        sofs[obji,6] = binskel.sum()
        sofs[obji,7] = hfeats[0]
        sofs[obji,8] = sofs[obji,6]/sofs[obji,0]
        sofs[obji,9] = (binobj*protimg).sum()/(binskel*protobj).sum()
        sofs[obji,10] = no_of_branch_points/sofs[obji,6]
    return sofs

