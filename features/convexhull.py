from __future__ import division
from matplotlib.nxutils import points_inside_poly, pnpoly
from ..imageprocessing.bbox import bbox
from numpy import *
import warnings

__all__ = ['convexhull','computeConvexHull']

def convexhull(bwimg):
    """
    hull = convexhull(bwimg)

    Compute the convex hull of the binary image and return it as a binary mask
    """
    # This function is still a bottleneck.
    # Computing the convex hull (computeConvexHull) is now very fast (in C),
    # but iterating over all possible points and calling pnpoly() is slow
    # a more specialised approach (like a line sweep) would do better
    Y,X=where(bwimg)
    P=list(zip(Y,X))
    if len(P) <= 3:
        return bwimg
    H=computeConvexHull(P)
    r,c=bwimg.shape
    min1,max1,min2,max2 = bbox(bwimg)
    res=zeros((r,c))
    for y in xrange(min1,max1+1):
        for x in xrange(min2,max2+1):
            res[y,x]=pnpoly(y,x,H)
    return res


def _isLeft(p0,p1,p2):
    """
    tests if a point is Left|On|Right of an infinite line.
    Input:  three points p0, p1, and p2

    Return:
        >0 for p2 left of the line through p0 and p1
        =0 for p2 on the line
        <0 for p2 right of the line
    """
    # Copyright 2001, softSurfer (www.softsurfer.com)
    # This code may be freely used and modified for any purpose
    # providing that this copyright notice is included with it.
    # SoftSurfer makes no warranty for this code, and cannot be held
    # liable for any real or imagined damage resulting from its use.
    # Users of this code must verify correctness for their application.
    return (p1[0]-p0[0])*(p2[1]-p0[1]) - (p2[0]-p0[0])*(p1[1]-p0[1])

def _inPlaceScan(P,reverse):
    P.sort(reverse=reverse)
    h=1
    N=len(P)
    for i in xrange(1,N):
        while h >= 2 and _isLeft(P[h-2],P[h-1],P[i]) >= 0:
            h -= 1
        t=P[i]
        P[i]=P[h]
        P[h]=t
        h += 1
    return h
def computeConvexHull(P):
    """
    From ``Space-Efficient Planar Convex Hull Algorithms''
        by Bronnimann et al.
    """
    try:
        import _convexhull
        return _convexhull.computeConvexHull(P)
    except Exception, e:
        warnings.warn('C code for convex hull failed. Resorting to (slow) python (Error: %s)' % e)
        h=_inPlaceScan(P,False)
        for i in xrange(h-1):
            t=P[i]
            P[i]=P[i+1]
            P[i+1]=t
        P2=P[h-2:]
        h_=_inPlaceScan(P2,True)
        return P[:h]+P2[:h_]

# vim: set ts=4 sts=4 sw=4 expandtab smartindent:
