# 10 Aug 98 - M.V. Boland
# 11 Jul 01 - Moment calculations optimized - G. Porreca
# 08 Aug 01 - Modified to call C implementation of previous
#             optimizations - G. Porreca
# Jun 2, 2002 - M.Velliste: added SLF names
# Feb 7 2008 - L-P Coelho: Ported to Python

# Copyright (C) 2006-2009 Murphy Lab
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
import numpy as np
from mahotas import euler
from scipy import ndimage

__all__ = ['imgfeatures','imgfeaturesdna']

def _norm2(V):
    return np.sqrt( np.dot(V,V) )

def imgfeatures(imageproc, isfield):
    """
    values = imgfeatures(imageproc, isfield=False)

    Calculates object features for imageproc
    where imageproc contains the pre-processed fluorescence image, 
    Pre-processed means that the image has been cropped and had 
    pixels of interest selected (via a threshold, for instance).
 
    Features calculated include:
      - Number of objects
      - Euler number of the image (# of objects - # of holes)
      - Average of the object sizes
      - Variance of the object sizes
      - Ratio of the largest object to the smallest
      - Average of the object distances from the COF [unless isfield]
      - Variance of the object distances from the COF [unless isfield]
      - Ratio of the largest object distance to the smallest

    @see imgfeaturesdna
    """

    return imgfeaturesdna(imageproc, None, isfield)

def imgfeaturesdna(imageproc, dnaproc, isfield=False):
    """
    values = imgfeaturesdna(imageproc, dnaproc, isfield=False)

    calculates object features for imageproc
    where imageproc contains the pre-processed fluorescence image, 
    and dnaproc the pre-processed dna fluorescence image.
    Pre-processed means that the image has been cropped and had 
    pixels of interest selected (via a threshold, for instance).
    use dnaproc=None to exclude features based on the dna image.  
    
    Parameters
    ----------
        * imageproc: pre-processed image
        * dnaproc: pre-processed DNA image
        * isfield: whether to calculate only field level features
                (default: False)

    Features calculated include:
      - Number of objects
      - Euler number of the image (# of objects - # of holes)
      - Average of the object sizes
      - Variance of the object sizes
      - Ratio of the largest object to the smallest
      - Average of the object distances from the COF [unless isfield]
      - Variance of the object distances from the COF [unless isfield]
      - Ratio of the largest object distance to the smallest
      - DNA: average of the object distances to the DNA COF [unless isfield]
      - DNA: variance of the object distances to the DNA COF [unless isfield]
      - DNA: ratio of the largest object distance to the smallest
      - DNA/Image: distance of the DNA COF to the image COF [unless isfield]
      - DNA/Image: ratio of the DNA image area to the image area
      - DNA/Image: fraction of image that overlaps with DNA 
    """
    bwimage = (imageproc > 0)
    imagelabeled,obj_number = ndimage.label(bwimage)
    values = [obj_number]
    dnacof = None

    if obj_number == 0:
        if isfield:
            if dnaproc is None: return np.zeros(5)
            raise NotImplementedError
        if dnaproc is None: return np.zeros(8)
        return np.zeros(14)

    euler_nr = euler(bwimage)
    values.append(euler_nr)

    cof = None
    if not isfield:
        cof = np.array(ndimage.center_of_mass(imageproc))

    # Find the maximum and minimum object sizes, and the distance 
    #    of each object to the center of fluorescence
    obj_distances = []
    obj_dnadistances = []

    moment_length = obj_number  + 1
    img_moment00 = np.zeros(moment_length);
    img_moment10 = np.zeros(moment_length);
    img_moment01 = np.zeros(moment_length);
    obj_sizes = np.zeros(moment_length);
    
    r,c=imageproc.shape
    try:
        from scipy import weave
        from scipy.weave import converters
        code='''
        for (int x = 0; x != c; ++x) {
            for (int y = 0; y != r; ++y) {
                unsigned moment_array_index = imagelabeled(y,x);
                img_moment00(moment_array_index) += imageproc(y,x);
                img_moment10(moment_array_index) += (x * imageproc(y,x));
                img_moment01(moment_array_index) += (y * imageproc(y,x));
                ++obj_sizes(moment_array_index);
            }
        }
        '''
        weave.inline(code,
            ['r','c','imageproc','imagelabeled','obj_sizes','img_moment00','img_moment01','img_moment10'],
            type_converters=converters.blitz)
    except ImportError:
        img_moment00 = np.zeros(moment_length);
        img_moment10 = np.zeros(moment_length);
        img_moment01 = np.zeros(moment_length);
        obj_sizes = np.zeros(moment_length);
        for x in xrange(c):
            for y in xrange(r):
                moment_array_index = imagelabeled[y,x]
                img_moment00[moment_array_index] += imageproc[y,x];
                img_moment10[moment_array_index] += (x * imageproc[y,x]);
                img_moment01[moment_array_index] += (y * imageproc[y,x]);
                obj_sizes[moment_array_index] += 1

    obj_sizes = np.array(obj_sizes[1:]) # Ignore object 0, i.e. background
    obj_size_avg = obj_sizes.mean()
    obj_size_var = obj_sizes.var()
    obj_size_ratio = obj_sizes.max()/obj_sizes.min() 

    values.append(obj_size_avg)
    values.append(obj_size_var) 
    values.append(obj_size_ratio)


    if not isfield:
        if dnaproc is not None:
            dnacof = np.array(ndimage.center_of_mass(dnaproc))
        for i in xrange(1,obj_number + 1):
            obj_m00 = float(img_moment00[i]);
            obj_m10 = float(img_moment10[i]);
            obj_m01 = float(img_moment01[i]);

            obj_center = np.array([obj_m01,obj_m10],float)/obj_m00;
            obj_distance = _norm2(obj_center - cof)
            
            obj_distances.append(obj_distance)

            if dnaproc is not None:
                obj_dnadistance = _norm2(obj_center - dnacof)
                obj_dnadistances.append(obj_dnadistance)
        obj_distances = np.array(obj_distances)
        obj_dist_avg = obj_distances.mean() ;
        obj_dist_var = obj_distances.var()
        mindist = obj_distances.min()
        if mindist != 0:
            obj_dist_ratio = obj_distances.max()/mindist
        else:
            obj_dist_ratio = 0

        values.append(obj_dist_avg)
        values.append(obj_dist_var)
        values.append(obj_dist_ratio)

        if dnaproc is not None:
            obj_dnadistances = np.array(obj_dnadistances)
            obj_dnadist_avg = obj_dnadistances.mean()
            obj_dnadist_var = obj_dnadistances.var()
            obj_dnamindist=obj_dnadistances.min()
            if obj_dnamindist != 0:
                obj_dnadist_ratio = obj_dnadistances.max()/obj_dnamindist 
            else:
                obj_dnadist_ratio = 0 ;
            values.append(obj_dnadist_avg)
            values.append(obj_dnadist_var)
            values.append(obj_dnadist_ratio)



    if dnaproc is not None:
        dna_area = (dnaproc > 0).sum()
        image_area = (imageproc > 0).sum()
        
        # what fraction of the image fluorescence area overlaps the dna image?
        image_overlap = ( (imageproc > 0) & (dnaproc > 0) ).sum()

        if image_area == 0:
            dna_image_area_ratio = 0
            image_dna_overlap = 0
        else:
            dna_image_area_ratio = dna_area/image_area ;
            image_dna_overlap = image_overlap/image_area ;
        
        if dnacof is not None:
            dna_image_distance = _norm2(cof-dnacof)
            values.append(dna_image_distance)
        values.append(dna_image_area_ratio)
        values.append(image_dna_overlap)
    return values

imgfeatures.names = [
    'object:number',
    'object:EulerNumber',
    'object_size:average',
    'object_size:variance',
    'object_size:ratio',
    'object_distance:average',
    'object_distance:variance',
    'object_distance:ratio']
imgfeatures.slfnames = [
    'SLF1.1',
    'SLF1.2',
    'SLF1.3',
    'SLF1.4',
    'SLF1.5',
    'SLF1.6',
    'SLF1.7',
    'SLF1.8']


imgfeaturesdna.names = imgfeatures.names + [
    'DNA_object_distance:average',
    'DNA_object_distance:variance',
    'DNA_object_distance:ratio',
    'DNA/image:distance',
    'DNA/image:area_ratio',
    'DNA/image:overlap']

imgfeaturesdna.slfnames = imgfeatures.slfnames + [
    'SLF2.17',
    'SLF2.18',
    'SLF2.19',
    'SLF2.20',
    'SLF2.21',
    'SLF2.22']
# vim: set ts=4 sts=4 sw=4 expandtab smartindent:
