from numpy import *
LEFT = 1
RIGHT = 2
UP = 4
DOWN = 8
ALL_DIRECTIONS = (LEFT|RIGHT|UP|DOWN)

def haralickfeatures(img,directions = ALL_DIRECTIONS):
    """
    Computes Haralick texture features for img.

    Returns a 13-vector of doubles
    """
    p=computecooccurence(img,directions)
    N=256
    feats=zeros(13)
    px=p.sum(0)
    py=p.sum(1)
    k=arange(N)
    ux=(k*px).sum()
    uy=(k*py).sum()
    sx=sqrt( (px*(k-ux)**2).sum() )
    sy=sqrt( (py*(k-uy)**2).sum() )
    px_plus_y=array([p.trace(G-N) for G in xrange(2*N)])
    px_minus_y=array([p.trace(G)+p.trace(-G) for G in xrange(N)])
    px_minus_y[0] /= 2
    g=array([[abs(i-j) for i in xrange(N)] for j in xrange(N)])
    i,i=mgrid[:N,:N]

    feats[0]=(p**2).sum()
    feats[1]=(k**2*px_minus_y).sum()
    feats[2]=1./sx/sy * ((i*j*p).sum() - ux*uy)
    feats[3]=((k-ux)**2*px).sum()
    feats[4]=(1./(1+(i-j)**2)*p).sum()
    feats[5]=(arange(2*N)*px_plus_y).sum()
    feats[6]=((arange(2*N)-feats[5])**2*px_plus_y).sum() # There is some confusion w.r.t. this feature.
                                                         # In some sources, it's feats[7] that is used
    feats[7]=entropy(px_plus_y)
    feats[8]=entropy(p)
    feats[9]=px_minus_y.var()
    feats[10]=entropy(px_minus_y)
    
    HX=entropy(px)
    HY=entropy(py)
    crosspxpy=dot(reshape(px,(N,1)),reshape(py,(1,N)))
    crosspxpy[crosspxpy == 0]=1. # This makes the log be zero and everything works OK below:
    HXY1=-(p*log(crosspxpy)).sum()
    HXY2=entropy(crosspxpy)

    feats[11]=feats[8]-HXY1/max(HX,HY)
    feats[12]=sqrt(1-exp(-2*(HXY2-feats[8])))

    return feats

def entropy(p):
    p=p[p != 0]
    return -(p*log2(p)).sum()

def computecooccurence(img,directions):
    comap = zeros((256,256))
    N,M=img.shape
    for i in xrange(N):
        for j in xrange(M):
            p=img[i,j]
            if directions & LEFT and i > 0:
                n=img[i-1,j]
                comap[p,n] += 1
            if directions & UP and j > 0:
                n=img[i,j-1]
                comap[p,n] += 1
            if directions & RIGHT and i < N -1:
                n=img[i+1,j]
                comap[p,n] += 1
            if directions & DOWN and j < M - 1:
                n=img[i,j+1]
                comap[p,n] += 1
    comap /= comap.sum()
    return comap


# vim: set ts=4 sts=4 sw=4 expandtab smartindent:
