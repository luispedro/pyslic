from matplotlib.nxutils import points_inside_poly
from numpy import *

__all__ = ['convexhull']

def convexhull(bwimg):
    """
    hull = convexhull(bwimg)

    Compute the convex hull of the binary image and return it as a binary mask
    """
    X,Y=where(bwimg)
    P=[(x,y) for x,y in zip(X,Y)]
    if len(P) < 2:
        return bwimg
    H=computesConvexHull(P)
    X,Y=bwimg.shape
    X,Y=mgrid[:X,:Y]
    res=points_inside_poly(array([(x,y) for x,y in zip(X.ravel(),Y.ravel())]),array(H))
    res=reshape(res,bwimg.shape)
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
def computesConvexHull(P):
    """
    From ``Space-Efficient Planar Convex Hull Algorithms''
        by Bronnimann et al.
    """
    h=_inPlaceScan(P,False)
    for i in xrange(h-1):
        t=P[i]
        P[i]=P[i+1]
        P[i+1]=t
    P2=P[h-2:]
    h_=_inPlaceScan(P2,True)
    return P[:h]+P2[:h_]

# vim: set ts=4 sts=4 sw=4 expandtab smartindent:
