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
    // isLeft(): tests if a point is Left|On|Right of an infinite line.
    //    Input:  three points P0, P1, and P2
    //    Return: >0 for P2 left of the line through P0 and P1
    //            =0 for P2 on the line
    //            <0 for P2 right of the line
    """
    # Copyright 2001, softSurfer (www.softsurfer.com)
    # This code may be freely used and modified for any purpose
    # providing that this copyright notice is included with it.
    # SoftSurfer makes no warranty for this code, and cannot be held
    # liable for any real or imagined damage resulting from its use.
    # Users of this code must verify correctness for their application.
     
    return (p1[0]-p0[0])*(p2[1]-p0[1]) - (p2[0]-p0[0])*(p1[1]-p0[1])

def computesConvexHull(P):
    """
    Compute convex hull based on Graham's scan algorithm as explained in
    http://www.softsurfer.com/Archive/algorithm_0109/algorithm_0109.htm
    """
    p0=P[0]
    for p in P:
        if p[0] < p0[0] or (p[0] == p0[0] and p[1] < p0[1]):
            p0=p
    P1=P[1:]
    def left(p1,p2):
        v=_isLeft(p0,p1,p2)
        if v > 0: return -1
        if v < 0: return 1
        return 0
    P1.sort(left) 
    p1=P1[0]
    stack=[p0,p1]
    i=1
    N=len(P1)
    while i < N:
        p1=stack[-1]
        p0=stack[-2]
        if _isLeft(p0,p1,P1[i]) > 0:
            stack.append(P1[i])
            i += 1
        else:
            del stack[-1]

    return stack

# vim: set ts=4 sts=4 sw=4 expandtab smartindent:
