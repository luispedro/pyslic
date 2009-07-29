'''
Created on Jul 15, 2009

@author: Robert Webb
'''
import numpy as np
def LBP(ImgArray, radius, point):
    '''
    LBP(image, radius, points)
    Imagelocation is the raw image to be evaluated.
    
    Firstly, the image is processed into an ndarray.
    The code first computes the points to be compared to the first pixel using bilinear interpolation
    Then, it compares these values to the center point, giving a 1 if higher or a 0 if lower.
    These values are put into an array.
    The array is rolled, and each value of the roll is put into a list.
    These values get turned into strings, which are then sorted
    This gives the lowest possible number created from interpolation
    Then the next pixel is checked.
    If it is the end of a row, the next row is checked.
    If it is the last pixel, the features are put into a dictionary with the number of times seen.
    '''
    angle = (2*np.pi)/point    
    feature = []
    lis = []
    ranging = range(2**(point))
    def binary(i,n):
        j =0
        return list((0,1)[i>>j & 1] for j in xrange(n-1,-1,-1))
    for i in ranging:
        q = ranging[i]
        lis = binary(q,q)
        if len(lis) < point:
            for i in xrange(point-len(lis)):
                lis.insert(0,0)
        if len(lis) > point:
            for i in xrange(len(lis)-point):
                del lis[0]
        feature.append(lis)
        lis = []
    final = [0] * (len(feature))
        
    for row in xrange(radius, ImgArray.shape[0]-radius-1):
        for pix in xrange(radius, ImgArray.shape[1]-radius-1):
            sign = [] # For the comparison function
            center = ImgArray[row,pix]
     
            for l in xrange(point): # For every point calculate these features
                x = radius * np.sin(angle*l)
                y = radius * np.cos(angle*l)
                ix = int(x)
                iy = int(y)
                a = ImgArray[(row + iy),(pix + ix)]
                b = ImgArray[(row + iy + 1),(pix + ix)]
                c = ImgArray[(row + iy),(pix + ix +1)]
                d = ImgArray[(row + iy+1),(pix + ix+1)]
                dx = x-ix
                dy = y-iy
                e = a+(c-a)*(dy)
                f = b+(d-b)*(dy)
                r = e+(f-e)*(dx)
                sign.append(int(center > r))
            sign = np.array(sign)
            bestval = sign.copy()

            for n in xrange(len(sign)):
                cur = np.roll(sign, n)
                for curbit, bestbit in zip(cur,bestval):
                    if curbit != bestbit:
                        if curbit < bestbit:
                            bestval = cur
                        break
            bestval = list(bestval)
            g = feature.index((bestval))
            final[g] = final[g] + 1
    return final