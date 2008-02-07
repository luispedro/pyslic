# Created by MV 1/18/02
# Modified MV 6/2/02: Added SLF names

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

from __future__ import division
from numpy import *
from scipy.ndimage import convolve, label
from bbox import croptobbox
from mmthin import mmthin
from convexhull import convexhull

__all__ = ['imgskelfeats']

def imgskelfeats(imageproc):
    """
[NAMES, VALUES, SLFNAMES] = ML_IMGSKELFEATURES(IMAGEPROC) calculates 
   skeleton features for IMAGEPROC

   where IMAGEPROC contains the pre-processed fluorescence image, 
   Pre-processed means that the image has been cropped and had 
   pixels of interest selected (via a threshold, for instance).
    """
    values = [] ;

    # Find objects in the image
    imagelabeled,N = label(imageproc > 0)
    for i in xrange(1,N+1):
        # Get an image of the single object
        objmask = (imagelabeled==i)
        objimage = imageproc * objmask
        # Compute skeleton features
        skelfeats=objskelfeats(objimage)
        values.append(skelfeats)

    # Average the skeleton features over the whole cell
    values=array(values)
    values=mean(values)
    # SLF names
    slfnames = []
    for feat_no in xrange(1,6):
        slfnames.append('SLF7.%s' % (feat_no+79))
    return values


def find_branch_points(img):
    """
    branch_points = ml_find_branch_points( img)
    img is the image with the skeleton of one object
    bran_points is the image with the branch points of the skeleton of that object
    """
    kernel = array([[0,1,0],[1,0,1],[0,1,0]])
    branch_points = convolve(img,kernel)*img;
    return (branch_points < 3) * (branch_points > 0)

def objskelfeats(objimg):
    """
[FEATS, NAMES] = ML_OBJSKELFEATS( OBJIMG)

Calculate skeleton features for the object OBJIMG.
    """
    objimg = croptobbox(objimg)
    objbin = objimg > 0
    objskel = mmthin(objbin);

    skellen = (objskel > 0).sum()
    objsize = objbin.sum()
    print objsize
    if objsize == 0:
        return array([0,0,0,0,0])

    skel_obj_area_ratio = skellen / objsize

    skelhull = convexhull(objskel);
    hullsize = skelhull.sum()

    # if hull size comes out smaller than length of skeleton then it
    # is obviously wrong, therefore adjust
    if hullsize < skellen:
         hullsize = skellen

    skel_hull_area_ratio = skellen / hullsize
    skel_fluor = objimg[objskel].sum()
    obj_fluor = objimg.sum()
    skel_obj_fluor_ratio = skel_fluor/obj_fluor

    branch_points = find_branch_points(objskel)
    no_of_branch_points = branch_points.sum()
    return array([skellen,skel_hull_area_ratio,skel_obj_area_ratio, skel_obj_fluor_ratio,no_of_branch_points/skellen])
    #    names = {'obj_skel_len' ...
    #            'obj_skel_hull_area_ratio' ...
    #            'obj_skel_obj_area_ratio' ...
    #            'obj_skel_obj_fluor_ratio' ...
    #            'obj_skel_branch_per_len'};
        


# vim: set ts=4 sts=4 sw=4 expandtab smartindent:
