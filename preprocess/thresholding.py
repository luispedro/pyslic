from numpy import *
from scipy.ndimage import histogram
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
    mint=0
    while hist[mint] == 0:
        mint += 1
        if mint == N:
            return mint 

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
    return res
        

def otsu(img):
    """
    T = otsu(img)

    Calculate a threshold according to the Otsu method.
    """
    raise Exception('Not implemented')

# vim: set ts=4 sts=4 sw=4 expandtab smartindent:
