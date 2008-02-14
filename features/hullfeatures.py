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

from __future__ import division
from bwperim import bwperim
from imgmoments import imgcentmoments
from convexhull import convexhull
from numpy import *

__all__ = ['hullfeatures']

def _bwarea(img):
    return (img > 0).sum()

def hullfeatures(imageproc,imagehull=None):
    """
    Compute hull features:

    hullfract:          bwarea/hullarea
    hullshape:          roundness of hull
    hull_eccentricity:  eccentricity of hull ellipse
    """

    if imagehull is None:
        imagehull = convexhull(imageproc > 0)

    Ahull = _bwarea(imagehull) ;
    hullfract = double((imageproc > 0).sum())/Ahull
    Phull = _bwarea(bwperim(imagehull))
    hullshape = (Phull**2)/(4*pi*Ahull)

    hull_mu00 = imgcentmoments(imagehull,0,0) ;
    hull_mu11 = imgcentmoments(imagehull,1,1) ;
    hull_mu02 = imgcentmoments(imagehull,0,2) ;
    hull_mu20 = imgcentmoments(imagehull,2,0) ;

# Parameters of the 'image ellipse'
#   (the constant intensity ellipse with the same mass and
#   second order moments as the original image.)
#   From Prokop, RJ, and Reeves, AP.  1992. CVGIP: Graphical
#   Models and Image Processing 54(5):438-460
    hull_semimajor = sqrt((2 * (hull_mu20 + hull_mu02 + \
                    sqrt((hull_mu20 - hull_mu02)**2 + \
                    4 * hull_mu11**2)))/hull_mu00) 

    hull_semiminor = sqrt((2 * (hull_mu20 + hull_mu02 - \
                    sqrt((hull_mu20 - hull_mu02)**2 + \
                    4 * hull_mu11**2)))/hull_mu00) 
    hull_eccentricity = sqrt(hull_semimajor**2 - hull_semiminor**2) / hull_semimajor ;
    #names = [cellstr('convex_hull:fraction_of_overlap') cellstr('convex_hull:shape_factor') cellstr('convex_hull:eccentricity')] ;
    #slfnames = [cellstr('SLF1.14') cellstr('SLF1.15') cellstr('SLF1.16')] ;
    #values = [hullfract hullshape hull_eccentricity] ;
    return array([hullfract,hullshape,hull_eccentricity])

# vim: set ts=4 sts=4 sw=4 expandtab smartindent:
