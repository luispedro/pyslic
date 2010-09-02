# Created by MV 1/18/02
# Modified MV 6/2/02: Added SLF names
# vim: set ts=4 sts=4 sw=4 expandtab smartindent:

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
import numpy
from scipy import ndimage
from mahotas.bbox import croptobbox
from mahotas import thin
from mahotas.polygon import fill_convexhull as convexhull

__all__ = ['imgskelfeatures', 'find_branch_points']

def imgskelfeatures(protproc):
    """
    values = imgskelfeatures(protproc)
    Compute skeleton features for protproc

    where protproc contains the pre-processed fluorescence image, 
    Pre-processed means that the image has been cropped and had 
    pixels of interest selected (via a threshold, for instance).
    """
    values = []

    # Find objects in the image
    labeled,N = ndimage.label(protproc > 0)
    if N == 0:
        return numpy.zeros(5,numpy.float64)
    objects = ndimage.find_objects(labeled)
    for i,slice in enumerate(objects):
        # Get an image of the single object
        objimage = protproc[slice] * (labeled[slice] == (i+1))
        # Compute skeleton features
        skelfeats=_objskelfeats(objimage)
        values.append(skelfeats)

    # Average the skeleton features over the whole cell
    values=numpy.array(values)
    values=values.mean(0)
    return values

imgskelfeatures.slfnames = [
    'SLF7.81',
    'SLF7.82',
    'SLF7.83',
    'SLF7.84',
    'SLF7.85',
    'SLF7.86']


def find_branch_points(img):
    """
    branch_points = find_branch_points( img)

    img is the image with the skeleton of one object
    branch_points is the image with the branch points of the skeleton of that object

    Ported from ml_find_branch_points.m
    """
    kernel = numpy.array([[0,1,0],[1,0,1],[0,1,0]])
    img=numpy.asarray(img,numpy.uint8)
    branch_points = img*ndimage.convolve(img,kernel,mode='constant')
    return (branch_points >= 3)

def _objskelfeats(objimg):
    """
    feats = _objskelfeats(objimg)

    Calculate skeleton features for the object OBJIMG.
    """
    objimg = objimg
    objbin = objimg > 0
    objsize = objbin.sum()

    if objsize == 0:
        return numpy.zeros(5)

    objskel = thin(objbin)
    skellen = objskel.sum()


    skelhull = convexhull(objskel);
    hullsize = skelhull.sum()
    hullsize = max(hullsize, skellen) # Corner cases such as [[1]]

    skel_hull_area_ratio = skellen / hullsize

    skel_obj_area_ratio = skellen/objsize

    skel_fluor = (objimg * objskel).sum()
    obj_fluor = objimg.sum()
    skel_obj_fluor_ratio = skel_fluor/obj_fluor

    branch_points = find_branch_points(objskel)
    no_of_branch_points = branch_points.sum()
    return numpy.array([skellen,skel_hull_area_ratio,skel_obj_area_ratio, skel_obj_fluor_ratio,no_of_branch_points/skellen])

imgskelfeatures.names = [
    'obj_skel_len',
    'obj_skel_hull_area_ratio',
    'obj_skel_obj_area_ratio',
    'obj_skel_obj_fluor_ratio',
    'obj_skel_branch_per_len']
        


