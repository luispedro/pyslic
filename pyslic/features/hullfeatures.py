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
import numpy
from mahotas import bwperim
from imgmoments import imgcentmoments
from mahotas.polygon import fill_convexhull as convexhull

from numpy import *
import numpy as np
from scipy.ndimage import center_of_mass

__all__ = ['hullfeatures','hullsizefeatures']

def _bwarea(img):
    if img.dtype != np.bool:
        img = (img > 0)
    return img.sum()

def _hull_computations(imageproc,imagehull = None):
    # Just share code between the two functions below
    if imagehull is None:
        imagehull = convexhull(imageproc > 0)

    Ahull = _bwarea(imagehull)
    Phull = _bwarea(bwperim(imagehull))

    cofy,cofx = center_of_mass(imagehull)
    hull_mu00 = imgcentmoments(imagehull,0,0,cofy,cofx)
    hull_mu11 = imgcentmoments(imagehull,1,1,cofy,cofx)
    hull_mu02 = imgcentmoments(imagehull,0,2,cofy,cofx)
    hull_mu20 = imgcentmoments(imagehull,2,0,cofy,cofx)

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
    return imagehull,Ahull, Phull, hull_semimajor, hull_semiminor

def hullsizefeatures(imageproc,imagehull=None):
    '''
    Compute image size features

    Area
    sqrt(Area)
    Perimeter
    Hull Semi-major axis
    Hull Semi-minor axis
    '''

    _,Ahull, Phull, hull_semimajor, hull_semiminor = _hull_computations(imageproc,imagehull)
    values=array([Ahull,sqrt(Ahull),Phull,hull_semimajor,hull_semiminor])
    values=r_[values,1./(values+(values==0))]
    return values

hullsizefeatures.names = [
'hullsize:area',
'hullsize:sqrt(area)',
'hullsize:perimeter',
'hullsize:semimajor',
'hullsize:semiminor',

'hullsize:inv(area)',
'hullsize:inv(sqrt(area))',
'hullsize:inv(perimeter)',
'hullsize:inv(semimajor)',
'hullsize:inv(semiminor))'
]

def hullfeatures(imageproc,imagehull=None):
    """
    Compute hull features:

    hullfract:          bwarea/hullarea
    hullshape:          roundness of hull
    hull_eccentricity:  eccentricity of hull ellipse
    """

    if imagehull is None:
        imagehull = convexhull(imageproc > 0)

    imagehull,Ahull, Phull, hull_semimajor, hull_semiminor = _hull_computations(imageproc,imagehull)
    if Ahull == 0: return numpy.array([0,0,0])
    hullfract = double((imageproc > 0).sum())/Ahull
    hullshape = (Phull**2)/(4*pi*Ahull)

    if hull_semimajor != 0.:
        hull_eccentricity = sqrt(hull_semimajor**2 - hull_semiminor**2) / hull_semimajor
    else:
        hull_eccentricity = 0.
    return numpy.array([hullfract,hullshape,hull_eccentricity])

hullfeatures.names = ['convex_hull:fraction_of_overlap', 'convex_hull:shape_factor', 'convex_hull:eccentricity']
hullfeatures.slf_names = ['SLF1.14', 'SLF1.15', 'SLF1.16']

# vim: set ts=4 sts=4 sw=4 expandtab smartindent:
