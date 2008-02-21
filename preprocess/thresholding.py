from numpy import *
from scipy.ndimage import histogram
from .basics import fullhistogram
__all__=['rc','otsu']
def rc(img,remove_zeros=False):
    """
    T = rc(img)
    
    Calculate a threshold according to the RC method.
    """
    mint=img.min()
    maxt=img.max()
    if maxt == mint:
        return maxt
    hist = histogram(img,mint,maxt,maxt-mint+1)
    if remove_zeros:
        hist[0]=0
    N=hist.size

    # Precompute most of what we need:
    sum1 = cumsum(arange(N) * hist)
    sum2 = cumsum(hist)
    sum3 = flipud(cumsum(flipud(arange(N) * hist)))
    sum4 = flipud(cumsum(flipud(hist)))

    maxt=N-1
    while hist[maxt] == 0:
        maxt -= 1

    res=maxt
    t=0
    while t < min(maxt,res):
        res=(sum1[t]/sum2[t] + sum3[t+1]/sum4[t+1])/2
        t += 1
    return res + mint
        

def otsu(img):
    """
    T = otsu(img)

    Calculate a threshold according to the Otsu method.
    """
    # Calculated according to http://homepages.inf.ed.ac.uk/rbf/CVonline/
    hist=fullhistogram(img)
    Ng=len(hist)
    nB=cumsum(hist)
    nO=nB[-1]-nB
    mu_B=0
    mu_O=(arange(1,Ng)*hist[1:]).sum()/hist.sum()
    best=nB[0]*nO[0]*(mu_B-mu_O)*(mu_B-mu_O)
    bestT=0

    for T in xrange(1,Ng):
        mu_B = (mu_B*nB[T-1] + T*hist[T-1]) / nB[T]
        mu_O = (mu_O*nO[T-1] - T*hist[T-1]) / nO[T]
        sigma_between=nB[T]*nO[T]*(mu_B-mu_O)*(mu_B-mu_O)
        if sigma_between > best:
            best = sigma_between
            bestT = T
    return bestT

# vim: set ts=4 sts=4 sw=4 expandtab smartindent:
