from __future__ import division
from numpy import *
from scipy.ndimage import distance_transform_edt,label, median_filter
from ..preprocess.thresholding import otsu
import nucleidetection

__all__ = ['voronoi','gvoronoi']
def voronoi(img,centers,distance='euclidean'):
    '''
    labeled = voronoi(img,centers,distance)

    distance can be one of
        'euclidean' : use euclidean distance (default)
        'manhatan'  : use manhatan distance
    '''
    labels=zeros_like(img)
    r,c=img.shape
    def seuclidean(y,x,cy,cx):
        dy=y-cy
        dx=x-cx
        return dy*dy+dx*dx
    def manhatan(y,x,cy,cx):
        dy=y-cy
        dx=x-cx
        return abs(dy)+abs(dx)
    dist=seuclidean
    if distance == 'manhatan' or distance == 'cityblock' or distance == 'city_block':
        dist=manhatan
    for y in xrange(r):
        for x in xrange(c):
            best=inf
            for i,(cy,cx) in enumerate(centers):
                now=dist(y,x,cy,cx)
                if now < best:
                    best=now
                    besti=i+1
            labels[y,x]=besti
    return labels

def gvoronoi(labeled,distance='euclidean'):
    """
    segmented = gvoronoi(labeled,labeled,distance)

    Generalised Voronoi Transform

    INPUT:
    * labeled: an array, of a form similar to the return of scipy.ndimage.label()
    * distance: one of 'euclidean', 'manhatan'

    RETURN
    segmented is of the same size and type as labeled and
        segmented[y,x] is the label of the object at position y,x
    """
    if distance == 'euclidean':
        L1,L2=distance_transform_edt(labeled== 0, return_distances=False,return_indices=True)
    else:
        raise Exception('gvoronoi: Distance "%s" not implemented' % distance)
    return labeled[L1,L2]

def gvoronoi_dna(dnaimg,distance='euclidean'):
    """
    labeled = gvoronoi_dna(dnaimg,distance)

    Generalised Voronoi Transform
    """
    labelednuclei,_=nucleidetection.labelnuclei(dnaimg)
    return gvoronoi(labelednuclei,distance)

# vim: set ts=4 sts=4 sw=4 expandtab smartindent:
