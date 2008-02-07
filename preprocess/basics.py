from scipy.ndimage import histogram
def fullhistogram(img):
    """
    Return a histogram containing all the values in the image
    """
    mint=img.min()
    maxt=img.max()
    return histogram(img,mint,maxt,maxt-mint+1)



# vim: set ts=4 sts=4 sw=4 expandtab smartindent:
