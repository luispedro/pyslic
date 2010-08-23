# -*- coding: utf-8 -*-
# Copyright (C) 2006--2008  Murphy Lab
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

# 19 Dec 98 - M.V. Boland
# Ported to Python by Luis Pedro Coelho <lpc@cmu.edu>

from __future__ import division
import mahotas.zernike

__all__ = ['zernike']

def zernike(img,D,radius,scale):
    """
    zvalues = zernike(img,D,radius) zernike moments through degree D

    Returns a vector of absolute Zernike moments through degree D for the
    image I.

    Parameters
    ----------
       * radius is used as the maximum radius for the Zernike polynomials.
       * scale is the scale of the image.

    Reference: Teague, MR. (1980). Image Analysis via the General
      Theory of Moments.  J. Opt. Soc. Am. 70(8):920-930.
    """
    return mahotas.zernike.zernike(img, D, radius/float(scale))


def znames(D,radius):
    """
     names = znames(D,radius) name oof zernike moments through degree D

    """
    names = []
    for n in xrange(D+1):
        for l in xrange(n+1):
            if (n-l)%2 == 0:
                names.append('Z_%s,%s' % (n, l))
    return names

# vim: set ts=4 sts=4 sw=4 expandtab smartindent:
