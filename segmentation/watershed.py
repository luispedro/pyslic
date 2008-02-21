from __future__ import division
from numpy import *
from heapq import *

__all__ = ['watershed']

def watershed(image,seeds):
    '''
    labeled = watershed(image,cofs)

    Basic, queue based, watershed algorithm
    '''
    queue=[]
    r,c=image.shape
    if type(seeds) == list:
        for i,(cofy, cofx) in enumerate(seeds):
            heappush(queue, (0,i+1,int(cofy),int(cofx)) )
    else:
        assert seeds.shape == (r,c)
        for y in xrange(r):
            for x in xrange(c):
                if seeds[y,x] > 0:
                    heappush(queue, (image[y,x],seeds[y,x],y,x))
    labeled=zeros_like(image)
    while queue:
        val,obj,y,x = heappop(queue)
        if labeled[y,x] != 0:
            continue
        labeled[y,x]=obj
        for i in xrange(3):
            ny=y-1+i
            if ny < 0 or ny >= r:
                continue
            for j in xrange(3):
                nx=x-1+j
                if nx < 0 or nx >= c:
                    continue
                heappush( queue, (image[ny,nx],obj,ny,nx) )
    return labeled

# vim: set ts=4 sts=4 sw=4 expandtab smartindent:
