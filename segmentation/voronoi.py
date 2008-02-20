from __future__ import division
from numpy import *
__all__ = ['voronoi']
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
    if distance == 'manhatan':
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
# vim: set ts=4 sts=4 sw=4 expandtab smartindent:
