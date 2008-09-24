# -*- coding: utf-8 -*-
# Copyright (C) 2008  Murphy Lab
# Carnegie Mellon University
# 
# Written by Luís Pedro Coelho <lpc@cmu.edu>
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published
# by the Free Software Foundation; either version 2 of the License,
# or (at your option) any later version.
#
# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
# 02110-1301, USA.
#
# For additional information visit http://murphylab.web.cmu.edu or
# send email to murphy@cmu.edu

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
