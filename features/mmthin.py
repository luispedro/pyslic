# Copyright (C) 2006  Murphy Lab
# Carnegie Mellon University
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

from numpy import *
from scipy.ndimage import binary_hit_or_miss
from scipy.misc.pilutil import imshow

__all__ = ['mmthin']
def mmthin(binimg):
    """
    skel = mmthin(binimg)
    image transformation by thinning

    rewrite the mmthin function in the morphological toolbox
    This code is written by yenixsa and Sam in Summer 2004
    Last updated on 12/3/2005
    Ported to Python by LPC on Feb 2008
    """

    binimg=binimg.copy()
    degrees = 45;
    num_elem = abs(360//degrees)

    struct_elem = []
    struct_elem.append([
            [0,0,0],
            [2,1,2],
            [1,1,1]])
    struct_elem.append([
            [2,0,0],
            [1,1,0],
            [1,1,2]])
    struct_elem.append([
            [1,2,0],
            [1,1,0],
            [1,2,0]])
    struct_elem.append([
            [1,1,2],
            [1,1,0],
            [2,0,0]])
    struct_elem.append([
            [1,1,1],
            [2,1,2],
            [0,0,0]])
    struct_elem.append([
            [2,1,1],
            [0,1,1],
            [0,0,2]])
    struct_elem.append([
            [0,2,1],
            [0,1,1],
            [0,2,1]])
    struct_elem.append([
            [0,0,2],
            [0,1,1],
            [2,1,1]])

    struct_elem=[array(E) for E in struct_elem]
    r,c=binimg.shape
    acnum_elem = 0;
    total_op=0

    while True:
        image_exp = zeros((r+2, c+2))
        image_exp[1:r+1, 1:c+1] = binimg
        newimg=hitmiss(image_exp,struct_elem[acnum_elem])
        binimg -= newimg[1:r+1,1:c+1]
        total_op += newimg.sum()

        acnum_elem +=  1;
        if acnum_elem == num_elem:
            if total_op == 0:
                break
            acnum_elem = 0
            total_op = 0
    return binimg

def hitmiss(binimg,struct_elem):
    '''
    Implementation of hit-or-miss operation
    '''
# Adapted from ml_mmhitmiss

    r,c=binimg.shape
    changed_image = zeros((r,c))
    try:
        from scipy import weave
        from scipy.weave import converters
        code = '''
#line 105 "mmthin.py"
        for (int y = 0; y != r-2 ; ++y) {
            for (int x = 0; x != c-2; ++x) {
                int hits = 0;
                for (int w = 0; w != 3; ++w) {
                    for (int z = 0; z != 3; ++z) {
                        if (struct_elem(w,z) == binimg(y+w,x+z) || struct_elem(w,z) == 2) ++hits;
                    }
                }
                if (hits == 9) changed_image(y+1,x+1) = 1;
            }
        }
        '''
        weave.inline(code,
            ['r','c','binimg','changed_image','struct_elem'],
            type_converters=converters.blitz
            )
    except Exception, e:
        print 'Weave failed. Resorting to (slow) python code'
        changed_image = zeros((r-2,c-2))
        for y in xrange(r-2):
            for x in xrange(c-2):
                hits = 0
                for w in xrange(3):
                    for z in xrange(3):
                        if struct_elem[w,z] == binimg[y+w,x+z] or struct_elem[w,z] == 2:
                            hits += 1
                if hits == 9:
                    changed_image[y+1,x+1]=1
    return changed_image
# vim: set ts=4 sts=4 sw=4 expandtab smartindent:
