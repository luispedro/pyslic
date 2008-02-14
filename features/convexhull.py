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

def _myDet(p, q, r):
    """Calc. determinant of a special matrix with three 2D points.

    The sign, "-" or "+", determines the side, right or left,
    respectivly, on which the point r lies, when measured against
    a directed vector from p to q.
    """

    # We use Sarrus' Rule to calculate the determinant.
    # (could also use the Numeric package...)
    sum1 = q[0]*r[1] + p[0]*q[1] + r[0]*p[1]
    sum2 = q[0]*p[1] + r[0]*q[1] + p[0]*r[1]

    return sum1 - sum2


def _isRightTurn((p, q, r)):
    "Do the vectors pq:qr form a right turn, or not?"
    assert p != q and q != r and p != r
    return _myDet(p, q, r) < 0

def computesConvexHull(P):
    """
	Calculate the convex hull of a set of points.
	Taken from http://aspn.activestate.com/ASPN/Cookbook/Python/Recipe/66527
    """

    # Get a local list copy of the points and sort them lexically.
    points=dict([(p,None) for p in P])
    points=points.keys()
    points.sort()

    # Build upper half of the hull.
    upper = [points[0], points[1]]
    for p in points[2:]:
	upper.append(p)
	while len(upper) > 2 and not _isRightTurn(upper[-3:]):
	    del upper[-2]

    # Build lower half of the hull.
    points.reverse()
    lower = [points[0], points[1]]
    for p in points[2:]:
	lower.append(p)
	while len(lower) > 2 and not _isRightTurn(lower[-3:]):
	    del lower[-2]

    # Remove duplicates.
    del lower[0]
    del lower[-1]

    # Concatenate both halfs and return.
    return tuple(upper + lower)



# vim: set ts=4 sts=4 sw=4 expandtab smartindent:
