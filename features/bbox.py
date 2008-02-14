from numpy import *
__all__ = ['bbox', 'croptobbox', 'expandto']
def bbox(img):
    """
    min1,max1,min2,max2 = bbox(img)

    Calculate the bounding box of image img.
    """
    min1=0
    min2=0
    max1,max2=img.shape
    while min1 < max1 and all(img[min1,:]==0):
        min1 += 1

    while max1 > min1 and all(img[max1-1,:]==0):
        max1 -= 1

    while min2 < max2 and all(img[:,min2]==0):
        min2 += 1

    while max2 > min2 and all(img[:,max2-1]==0):
        max2 -= 1

    return min1,max1,min2,max2

def croptobbox(img):
    """
    Returns a version of img cropped to the image's bounding box
    """
    
    min1,max1,min2,max2 = bbox(img)
    return img[min1:max1,min2:max2]

def expandto(img,S1,S2):
    """
    Returns a zero padded version of img
    so that it is of size (S1,S2).

    Img will be on the top left of the result.
    """

    result=zeros((S1,S2),img.dtype)
    M,N=img.shape
    assert M <= S1
    assert N <= S2
    st1=(S1-M)//2
    st2=(S2-N)//2
    result[st1:(st1+M),st2:(st2+N)]=img
    return result

# vim: set ts=4 sts=4 sw=4 expandtab smartindent:
