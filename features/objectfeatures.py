from __future__ import division
from ..image import Image
from numpy import *
from ..imageprocessing.bweuler import bweuler
from scipy.ndimage import *
from .hullfeatures import hullfeatures
from .mmthin import mmthin
from ..imageprocessing.convexhull import convexhull
from .imgskelfeats import find_branch_points
import warnings

def objectfeatures(protimg,dnaimg=None,objects=None):
    '''
    values=objectfeatures(img,labeled=None,objects=None)

    protimg is the preprocessed protein image
    dnaimg is the preprocessed dna image
        if not given, then dna dependent fields will be zero
    
    This implements the object features described in
    "Object Type Recognition for Automated Analysis of Protein Subcellular Location"
    by Ting Zhao, Meel Velliste, Michael V. Boland, and Robert F. Murphy
    in IEEE Transaction on Image Processing
    '''

    if type(protimg) is Image:
        if dnaimg is None:
            dnaimg = protimg.channeldata.get(Image.procdna_channel,None)
        else:
            warnings.warn('Bizarre combination of arguments for objectfeatures: both an Image and a separate DNA channel')
        protimg = protimg.channeldata[Image.procprotein_channel]

    labeled,N = label(protimg,ones((3,3)))
    if objects is None:
        objects = xrange(1,N+1)
    if type(objects) is int:
        objects = [objects]

    sofs = zeros((len(objects),11))
    if dnaimg is not None:
        dnacofy,dnacofx=center_of_mass(dnaimg)
        bindna=(dnaimg > 0)
    for obji,obj in enumerate(objects):
        binobj = (labeled == obj)
        objimg = protimg * binobj
        cofy,cofx = center_of_mass(objimg)
        objskel = mmthin(binobj);
        binskel = (objskel > 0)
        objhull=convexhull(binobj)
        no_of_branch_points = find_branch_points(objskel).sum()
        hfeats=hullfeatures(objimg,objhull)

        sofs[obji,0] = binobj.sum()
        if dnaimg is not None:
            sofs[obji,1] = sqrt((cofy-dnacofy)**2+(cofx-dnacofx)**2)
            sofs[obji,2] = (binobj*bindna).sum()/sofs[obji,0]
        sofs[obji,3] = hfeats[2]
        sofs[obji,4] = bweuler(binobj)
        sofs[obji,5] = hfeats[1]
        sofs[obji,6] = binskel.sum()
        sofs[obji,7] = hfeats[0]
        sofs[obji,8] = sofs[obji,6]/sofs[obji,0]
        sofs[obji,9] = (binobj*protimg).sum()/(binskel*protimg).sum()
        sofs[obji,10] = no_of_branch_points/sofs[obji,6]
    return sofs

