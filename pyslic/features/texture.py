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

__all__ = ['haralickfeatures','computecooccurence']

def haralickfeatures(img):
    """
    Computes Haralick texture features for img.
        

    Returns a 13-vector of doubles
    """
    if img.max() > 256:
        img = asarray(double(img) * 256./img.max(),uint8)

    nr_dirs = 4
    if len(img.shape) == 3:
        nr_dirs = 13

    feats=zeros((nr_dirs,13))
    for dir in xrange(nr_dirs):
        p=computecooccurence(img,dir)
        if p.size == 0: continue
        N,_=p.shape
        px=p.sum(0)
        py=p.sum(1)
        k=arange(N)
        ux=(k*px).sum()
        uy=(k*py).sum()
        sx=sqrt( (px*(k-ux)**2).sum() )
        sy=sqrt( (py*(k-uy)**2).sum() )
        px_plus_y=zeros(2*N)
        try:
            from scipy import weave
            from scipy.weave import converters
            px_minus_y = zeros(N)
            code = '''
            for (int i = 0; i != N; ++i) {
                for (int j = 0; j != N; ++j) {
                    px_plus_y(i+j) += p(i,j);
                    px_minus_y(std::abs(i-j)) += p(i,j);
                }
            }
            '''
            weave.inline(
                    code,
                    ['N','p','px_plus_y','px_minus_y'],
                    type_converters=converters.blitz)
        except:
            for i in xrange(N):
                for j in xrange(N):
                    px_plus_y[i+j] += p[i,j]
            px_minus_y=array([p.trace(G)+p.trace(-G) for G in xrange(N)])
            px_minus_y[0] /= 2

        i,j=mgrid[:N,:N]

        feats[dir,0]=(p**2).sum()
        feats[dir,1]=(k**2*px_minus_y).sum()

        feats[dir,2]=1./sx/sy * ((i*j*p).sum() - ux*uy)

        feats[dir,3]=((k-ux)**2*px).sum()
        feats[dir,4]=(1./(1+(i-j)**2)*p).sum()
        feats[dir,5]=(arange(2*N)*px_plus_y).sum()

        feats[dir,7]=entropy(px_plus_y)
        # There is some confusion w.r.t. this feature.
        # This is the formula in Haralick's paper, but some sources
        #  consider it a typo that feats[dir,7] is used and argue that the
        #  intended feature is obtained by substituting feats[dir,5].
        #
        # feats[dir,5] is probably the right one to use.
        # This version is the one in MurphyLab's Matlab (arguably incorrect) implementation
        feats[dir,6]=((arange(2*N)-feats[dir,7])**2*px_plus_y).sum()
                                                                 
        feats[dir,8]=entropy(p.ravel())
        feats[dir,9]=px_minus_y.var() # This is wrongly implemented in ml_texture
        feats[dir,10]=entropy(px_minus_y)
        
        HX=entropy(px)
        HY=entropy(py)
        crosspxpy=outer(px,py)
        crosspxpy[crosspxpy == 0]=1. # This makes the log be zero and everything works OK below:
        HXY1=-(p*log2(crosspxpy)).sum()
        HXY2=entropy(crosspxpy.ravel())

        feats[dir,11]=(feats[dir,8]-HXY1)/max(HX,HY)
        feats[dir,12]=sqrt(1-exp(-2*(HXY2-feats[dir,8])))

    return feats

haralickfeatures.names=[
    'angular_second_moment',
    'contrast','correlation',
    'sum_of_squares',
    'inverse_diff_moment',
    'sum_avg',
    'sum_var',
    'sum_entropy',
    'entropy',
    'diff_var',
    'diff_entropy',
    'info_measure_corr_1',
    'info_measure_corr_2']
haralickfeatures.slfnames = [ ('SLF3.%d' % n) for n in xrange(66,66+13)]

def entropy(p):
    try:
        from scipy import weave
        from scipy.weave import converters
        res=array([0.],double)
        p=p.ravel()
        N=len(p)
        code = '''
        for (int i = 0; i != N; ++i)
            if (p(i) > 0.) res(0) += - p(i) * log2(p(i));
        '''
        weave.inline(
                code,
                ['N','p','res'],
                type_converters=converters.blitz)
        return res[0]
    except:
        import scipy.stats
        p=p[p != 0]
        return scipy.stats.entropy(p)/log(2) # scipy.stats.entropy is natural log based!

_dir_deltas = [
    (1,0,0),
    (1,1,0),
    (0,1,0),
    (1,-1,0),
    (0,0,1),
    (1,0,1),
    (0,1,1),
    (1,1,1),
    (1,-1,1),
    (1,0,-1),
    (0,1,-1),
    (1,1,-1),
    (1,-1,-1),
]
def computecooccurence(img,dir,remove_zeros=True):
    assert dir >= 0
    assert dir < len(_dir_deltas)
    assert dir < 5 or len(img.shape) == 3
    if len(img.shape) == 2:
        img = img.reshape((1,)+img.shape)
    if img.min() < 0:
        img = img - img.min()
    Ng=img.max()+1
    comap = zeros((Ng,Ng))
    N0,N1,N2=img.shape
    dx,dy,dz=_dir_deltas[dir]
    try:
        from scipy import weave
        from scipy.weave import converters
        code = '''
#line 180 "texture.py"
        int nz,ny,nx;
        for (int z = 0; z != N0; ++z) {
            for (int y = 0; y != N1; ++y) {
                for (int x = 0; x != N2; ++x) {
                    nz=z+dz;
                    ny=y+dy;
                    nx=x+dx;
                    if (nz < 0 || nx < 0 || ny < 0 || nz >= N0 || ny >= N1 || nx >= N2) continue;
                    int p=img(z,y,x);
                    int n=img(nz,ny,nx);
                    ++comap(p,n);
                }
            }
        }
        '''
        weave.inline(code,
            ['N0','N1','N2','dz','dy','dx','comap','img'],
            type_converters=converters.blitz)
    except:
        for z in xrange(N0):
            for y in xrange(N1):
                for x in xrange(N2):
                    nz=z+dz
                    ny=y+dy
                    nx=x+dx
                    if nz < 0 or ny < 0 or nx < 0 or nz > N0 or ny >= N1 or nx >= N2:
                        continue
                    p=img[z,y,x]
                    n=img[nz,ny,nx]
                    comap[p,n] += 1
    comap = comap + comap.T
    if remove_zeros:
        comap=comap[1:,1:]
    comap /= comap.sum()
    return comap


# vim: set ts=4 sts=4 sw=4 expandtab smartindent:
