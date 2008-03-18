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
    try:
        from scipy import weave
        from scipy.weave import converters
        vals=array([max1,max2,min1,min2])
        code='''
        int max1=vals(0);
        int max2=vals(1);
        for (int y = 0; y != max1; ++y) {
            for (int x = 0; x != max2; ++x) {
                if (img(y,x) > 0) {
                    if (y < vals(0)) vals(0) = y;
                    if (x < vals(1)) vals(1) = x;
                    if (y > vals(2)) vals(2) = y;
                    if (x > vals(3)) vals(3) = x;
                }
            }
        }
        '''
        weave.inline(
                code,
                ['vals','img'],
                type_converters=converters.blitz)
        min1,min2,max1,max2=tuple(vals)
    except:
        pos1,=where(img.any(1))
        if len(pos1) == 0:
            return 0,0,0,0
        min1=pos1[0]
        max1=pos1[-1]

        pos2,=where(img.any(0))
        min2=pos2[0]
        max2=pos2[-1]
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
