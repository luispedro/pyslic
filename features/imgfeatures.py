# 10 Aug 98 - M.V. Boland
# 11 Jul 01 - Moment calculations optimized - G. Porreca
# 08 Aug 01 - Modified to call C implementation of previous
#             optimizations - G. Porreca
# Jun 2, 2002 - M.Velliste: added SLF names
# Feb 7 2008 - L-P Coelho: Ported to Python

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
from bweuler import bweuler
from scipy.ndimage import *

__all__ = ['imgfeatures']

def _norm2(V):
    return sqrt( (V**2).sum() )

def imgfeatures(imageproc,dnaproc):
    """
    values = imgfeatures(imageproc,dnaproc)

   Calculates features for IMAGEPROC
   where IMAGEPROC contains the pre-processed fluorescence image, 
   and DNAPROC the pre-processed DNA fluorescence image.
   Pre-processed means that the image has been cropped and had 
   pixels of interest selected (via a threshold, for instance).
   Use DNAPROC=[] to exclude features based on the DNA image.  

   Features calculated include:
     - Number of objects
     - Euler number of the image (# of objects - # of holes)
     - Average of the object sizes
     - Variance of the object sizes
     - Ratio of the largest object to the smallest
     - Average of the object distances from the COF
     - Variance of the object distances from the COF
     - Ratio of the largest object distance to the smallest
     - DNA: average of the object distances to the DNA COF
     - DNA: variance of the object distances to the DNA COF
     - DNA: ratio of the largest object distance to the smallest
     - DNA/Image: distance of the DNA COF to the image COF
     - DNA/Image: ratio of the DNA image area to the image area
     - DNA/Image: fraction of image that overlaps with DNA 
    """


    names = []
    slfnames = []

    bwimage=(imageproc > 0)
    imagelabeled,obj_number = label(bwimage)
    values = [obj_number]

    euler_nr = bweuler(bwimage)
    values.append(euler_nr)

    #names = [names cellstr('object:number') cellstr('object:EulerNumber')] ;
    #slfnames = [slfnames cellstr('SLF1.1') cellstr('SLF1.2')] ;

    # Calculate the center of fluorescence of IMAGE
    cof = array(center_of_mass(imageproc))

    dnacof=None
    if dnaproc is not None:
        dnacof=array(center_of_mass(dnaproc))
    
    # Find the maximum and minimum object sizes, and the distance 
    #    of each object to the center of fluorescence
    obj_distances = [] ;
    obj_dnadistances = [] ;


    ###################################################################################################
    # GP 08/08/01
    # Calls gp_imgmoments, which calls a C implementation of the code commented
    # out below; the C version runs approx. 100X faster (10,000#)
    #
    #[img_moment00, img_moment10, img_moment01, obj_size] = ml_moments_1(int32(imageproc),... 
    #                                     int32(imagelabeled));
    #
    #obj_sizes = obj_size(1,2:end);
    # /GP
    ###################################################################################################

    moment_length = obj_number  + 1
    img_moment00 = zeros(moment_length);
    img_moment10 = zeros(moment_length);
    img_moment01 = zeros(moment_length);
    obj_sizes = zeros(moment_length);
    
    
    r,c=imageproc.shape
    for x in xrange(c):
        for y in xrange(r):
            moment_array_index = imagelabeled[y,x]
            img_moment00[moment_array_index] += imageproc[y,x];
            img_moment10[moment_array_index] += (x * imageproc[y,x]);
            img_moment01[moment_array_index] += (y * imageproc[y,x]);
            obj_sizes[moment_array_index] += 1
    


    for i in xrange(1,obj_number + 1):
        obj_m00 = double(img_moment00[i]);
        obj_m10 = double(img_moment10[i]);
        obj_m01 = double(img_moment01[i]);

        obj_center = array([obj_m01,obj_m10],double)/obj_m00;
        obj_distance = _norm2(obj_center - cof)
        
        obj_distances.append(obj_distance)

        if dnaproc is not None:
            obj_dnadistance = _norm2(obj_center - dnacof)
            obj_dnadistances.append(obj_dnadistance)

    obj_sizes=array(obj_sizes[1:]) # Ignore object 0, i.e. background
    obj_size_avg = obj_sizes.mean()
    obj_size_var = obj_sizes.var()
    obj_size_ratio = obj_sizes.max()/obj_sizes.min() 

    #names = [names cellstr('object_size:average') ...
    #        cellstr('object_size:variance') ...
    #        cellstr('object_size:ratio')] ;
    #slfnames = [slfnames cellstr('SLF1.3') ...
    #        cellstr('SLF1.4') ...
    #        cellstr('SLF1.5')] ;
    values.append(obj_size_avg)
    values.append(obj_size_var) 
    values.append(obj_size_ratio)

    obj_distances=array(obj_distances)
    obj_dist_avg = mean(obj_distances) ;
    obj_dist_var = var(obj_distances)
    mindist=obj_distances.min()
    if mindist != 0:
        obj_dist_ratio = obj_distances.max()/mindist
    else:
        obj_dist_ratio = 0

    #names = [names cellstr('object_distance:average') ... 
    #        cellstr('object_distance:variance') ...
    #        cellstr('object_distance:ratio')] ;
    #slfnames = [slfnames cellstr('SLF1.6') ... 
    #        cellstr('SLF1.7') ...
    #        cellstr('SLF1.8')] ;
    values.append(obj_dist_avg)
    values.append(obj_dist_var)
    values.append(obj_dist_ratio)

    if dnaproc is not None:
        obj_dnadistances=array(obj_dnadistances)
        obj_dnadist_avg = obj_dnadistances.mean()
        obj_dnadist_var = obj_dnadistances.var()
        obj_dnamindist=obj_dnadistances.min()
        if obj_dnamindist != 0:
            obj_dnadist_ratio = obj_dnadistances.max()/obj_dnamindist 
        else:
            obj_dnadist_ratio = 0 ;
        #names = [names cellstr('DNA_object_distance:average') ...
        #                   cellstr('DNA_object_distance:variance') ...
        #                   cellstr('DNA_object_distance:ratio')] ;
        #slfnames = [slfnames cellstr('SLF2.17') ...
        #                   cellstr('SLF2.18') ...
        #                   cellstr('SLF2.19')] ;
        values.append(obj_dnadist_avg)
        values.append(obj_dnadist_var)
        values.append(obj_dnadist_ratio)

        dna_image_distance = _norm2(cof-dnacof)

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
        
        #names = [names cellstr('DNA/image:distance') ...
        #        cellstr('DNA/image:area_ratio') ...
        #        cellstr('DNA/image:overlap')] ;
        #slfnames = [slfnames cellstr('SLF2.20') ...
        #        cellstr('SLF2.21') ...
        #        cellstr('SLF2.22')] ;
        values.append(dna_image_distance)
        values.append(dna_image_area_ratio)
        values.append(image_dna_overlap)
    return values

# vim: set ts=4 sts=4 sw=4 expandtab smartindent:
