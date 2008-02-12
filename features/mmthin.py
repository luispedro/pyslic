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

__all__ = ['mmthin']
def mmthin(binimg):
    """
    ML_MMTHIN image transformation by thinning
    IMG_SKEL=ML_MMTHIN(BIN_IMAGE)

    rewrite the mmthin function in the morphological toolbox
    This code is written by yenixsa and Sam in Summer 2004
    Last updated on 12/3/2005
    Ported to Python by LPC on Feb 2008
    """

    binimg=binimg.copy()
    degrees = 45;
    num_elem = abs(360//degrees)

    struct_elem = []
    struct_elem.append([[0,0,0],[2,1,2],[1,1,1]])
    struct_elem.append([[2,0,0],[1,1,0],[1,1,2]])
    struct_elem.append([[1,2,0],[1,1,0],[1,2,0]])
    struct_elem.append([[1,1,2],[1,1,0],[2,0,0]])
    struct_elem.append([[1,1,1],[2,1,2],[0,0,0]])
    struct_elem.append([[2,1,1],[0,1,1],[0,0,2]])
    struct_elem.append([[0,2,1],[0,1,1],[0,2,1]])
    struct_elem.append([[0,0,2],[0,1,1],[2,1,1]])

    r,c=binimg.shape
    total_oper = 0;
    acnum_elem = 0;

    while True:
        acnum_elem +=  1;
        image_exp = zeros((r+4, c+4))
        image_exp[2:r+2, 2:c+2] = binimg
        newimg=binary_hit_or_miss(image_exp,struct_elem[acnum_elem])
        binimg -= newimg[2:r+2,2:c+2]

        if acnum_elem == num_elem - 1:
            acnum_elem = 0
            if newimg.sum() == 0:
                break
    return binimg

# vim: set ts=4 sts=4 sw=4 expandtab smartindent:
