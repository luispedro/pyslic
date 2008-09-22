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
from ..imageprocessing.bbox import croptobbox
from mmthin import mmthin
from ..imageprocessing.convexhull import convexhull

__all__ = ['imgskelfeatures', 'find_branch_points']

def imgskelfeatures(protproc):
    """
    values = imgskelfeatures(protproc)
    Compute skeleton features for protproc

   where protproc contains the pre-processed fluorescence image, 
   Pre-processed means that the image has been cropped and had 
   pixels of interest selected (via a threshold, for instance).
    """
    values = [] ;

    # Find objects in the image
    imagelabeled,N = label(protproc > 0)
    if N == 0:
        return zeros(5,double)
    for i in xrange(1,N+1):
        # Get an image of the single object
        objimage = protproc * (imagelabeled == i)
        # Compute skeleton features
        skelfeats=objskelfeats(objimage)
        values.append(skelfeats)

    # Average the skeleton features over the whole cell
    values=array(values)
    values=values.mean(0)
    # SLF names
    slfnames = []
    for feat_no in xrange(1,6):
        slfnames.append('SLF7.%s' % (feat_no+79))
    return values


def find_branch_points(img):
    """
    branch_points = find_branch_points( img)

    img is the image with the skeleton of one object
    branch_points is the image with the branch points of the skeleton of that object

    Ported from ml_find_branch_points.m
    """
    kernel = array([[0,1,0],[1,0,1],[0,1,0]])
    img=asarray(img,uint8)
    branch_points = img*convolve(img,kernel,mode='constant')
    return (branch_points >= 3)

def objskelfeats(objimg):
    """
    feats = objskelfeats(objimg)

    Calculate skeleton features for the object OBJIMG.
    """
    objimg = croptobbox(objimg)
    objbin = objimg > 0
    objsize = objbin.sum()

    if objsize == 0:
        return array([0,0,0,0,0])

    objskel = mmthin(objbin);
    skellen = (objskel > 0).sum()


    skelhull = convexhull(objskel);
    hullsize = skelhull.sum()

    # if hull size comes out smaller than length of skeleton then it
    # is obviously wrong, therefore adjust
    if hullsize < skellen:
         hullsize = skellen

    skel_hull_area_ratio = skellen / hullsize

    skel_obj_area_ratio = skellen/objsize

    skel_fluor = objimg[objskel].sum()
    obj_fluor = objimg.sum()
    skel_obj_fluor_ratio = skel_fluor/obj_fluor

    branch_points = find_branch_points(objskel)
    no_of_branch_points = branch_points.sum()
    return array([skellen,skel_hull_area_ratio,skel_obj_area_ratio, skel_obj_fluor_ratio,no_of_branch_points/skellen])

imgskelfeatures.names = [
    'obj_skel_len',
    'obj_skel_hull_area_ratio',
    'obj_skel_obj_area_ratio',
    'obj_skel_obj_fluor_ratio',
    'obj_skel_branch_per_len']
        


# vim: set ts=4 sts=4 sw=4 expandtab smartindent:
