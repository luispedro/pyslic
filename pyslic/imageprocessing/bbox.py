# -*- coding: utf-8 -*-
# Copyright (C) 2008  Murphy Lab
# vim: set ts=4 sts=4 sw=4 expandtab smartindent:
# Carnegie Mellon University
# 
# Written by Luis Pedro Coelho <lpc@cmu.edu>
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

from mahotas.bbox import bbox, croptobbox
import numpy as np
__all__ = ['bbox', 'croptobbox', 'expandto']

def expandto(img,S1,S2):
    """
    Returns a zero padded version of img
    so that it is of size (S1,S2).

    Img will be on the top left of the result.
    """

    result = np.zeros((S1,S2),img.dtype)
    M,N=img.shape
    assert M <= S1
    assert N <= S2
    st1=(S1-M)//2
    st2=(S2-N)//2
    result[st1:(st1+M),st2:(st2+N)]=img
    return result

