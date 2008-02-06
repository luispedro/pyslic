from numpy import *

__all__ = ['imgmoments','imgcentmoments']
def imgmonents(img,x,y):
    """
     M = imgmonents(IMG, X, Y) calculates the moment MXY for IMAGE
     imgmoments(img, x, y), 
        where IMG is the image to be processed and X and Y define
        the order of the moment to be calculated. For example, 

        imgmoments(image,0,1) calculates the first order moment 
        in the y-direction, and 

        imgmoments(image,0,1)/imgmoments(image,0,0) is the 
        'center of mass (fluorescence)' in the y-direction
    """
    img = double(img)
    N,M=img.shape
    img *= arange(M)**y
    img=img.T
    img *= arange(N)**x
    return img.sum()

def imgcentmoments(img,x,y):
    """
    M_xy = imgcentmoments(img,x,y)
    """
    img = double(img)
    cofx, cofy = center_of_mass(img)
    N,M=img.shape
    img *= (arange(M)-cofy)**y
    img = img.T
    img *= (arange(N)-cofx)**x
    return img.sum()

# vim: set ts=4 sts=4 sw=4 expandtab smartindent:
