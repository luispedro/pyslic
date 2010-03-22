# -*- coding: utf-8 -*-
# Copyright (C) 2008-2009  Murphy Lab
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

from __future__ import division, with_statement
from tas import pftasinfo
from overlap import overlapinfo

def _harprops(level):
    level -= 1
    START_DHAR = 2*9*2 + 1
    return [
        ('SLF33.%s' % (START_DHAR + level*13 +  0),'angular_second_moment-after-dsample-%s' % (level+1),2,1),
        ('SLF33.%s' % (START_DHAR + level*13 +  1),'contrast-after-dsample-%s' % (level + 1),2,1),
        ('SLF33.%s' % (START_DHAR + level*13 +  2),'correlation-after-dsample-%s' % (level + 1),2,1),
        ('SLF33.%s' % (START_DHAR + level*13 +  3),'sum_of_squares-after-dsample-%s' % (level + 1),2,1),
        ('SLF33.%s' % (START_DHAR + level*13 +  4),'inverse_diff_moment-after-dsample-%s' % (level + 1),2,1),
        ('SLF33.%s' % (START_DHAR + level*13 +  5),'sum_avg-after-dsample-%s' % (level + 1),2,1),
        ('SLF33.%s' % (START_DHAR + level*13 +  6),'sum_var-after-dsample-%s' % (level + 1),2,1),
        ('SLF33.%s' % (START_DHAR + level*13 +  7),'sum_entropy-after-dsample-%s' % (level + 1),2,1),
        ('SLF33.%s' % (START_DHAR + level*13 +  8),'entropy-after-dsample-%s' % (level + 1),2,1),
        ('SLF33.%s' % (START_DHAR + level*13 +  9),'diff_var-after-dsample-%s' % (level + 1),2,1),
        ('SLF33.%s' % (START_DHAR + level*13 + 10),'diff_entropy-after-dsample-%s' % (level + 1),2,1),
        ('SLF33.%s' % (START_DHAR + level*13 + 11),'info_measure_corr_1-after-dsample-%s' % (level + 1),2,1),
        ('SLF33.%s' % (START_DHAR + level*13 + 12),'info_measure_corr_2-after-dsample-%s' % (level + 1),2,1),
    ]

featinfo = [
    ('SLF7.80','obj_skel_len',2,1),
    ('SLF7.81','obj_skel_hull_area_ratio',2,1),
    ('SLF7.82','obj_skel_obj_area_ratio',2,1),
    ('SLF7.83','obj_skel_obj_fluor_ratio',2,1),
    ('SLF7.84','obj_skel_branch_per_len',2,1),
    ('SLF7.79','fract_nonobj_fluor',2,1),
    ('SLF1.1','object:number',2,1),
    ('SLF1.2','object:EulerNumber',2,1),
    ('SLF1.3','object_size:average',2,1),
    ('SLF1.4','object_size:variance',2,1),
    ('SLF1.5','object_size:ratio',2,1),
    ('SLF1.6','object_distance:average',2,1),
    ('SLF1.7','object_distance:variance',2,1),
    ('SLF1.8','object_distance:ratio',2,1),
    ('SLF2.17','DNA_object_distance:average',2,2),
    ('SLF2.18','DNA_object_distance:variance',2,2),
    ('SLF2.19','DNA_object_distance:ratio',2,2),
    ('SLF2.20','DNA/image:distance',2,2),
    ('SLF2.21','DNA/image:area_ratio',2,2),
    ('SLF2.22','DNA/image:overlap',2,2),
    ('SLF1.14','convex_hull:fraction_of_overlap',2,1),
    ('SLF1.15','convex_hull:shape_factor',2,1),
    ('SLF1.16','convex_hull:eccentricity',2,1),
    ('SLF3.17','Z_0_0',2,1),
    ('SLF3.18','Z_1_1',2,1),
    ('SLF3.19','Z_2_0',2,1),
    ('SLF3.20','Z_2_2',2,1),
    ('SLF3.21','Z_3_1',2,1),
    ('SLF3.22','Z_3_3',2,1),
    ('SLF3.23','Z_4_0',2,1),
    ('SLF3.24','Z_4_2',2,1),
    ('SLF3.25','Z_4_4',2,1),
    ('SLF3.26','Z_5_1',2,1),
    ('SLF3.27','Z_5_3',2,1),
    ('SLF3.28','Z_5_5',2,1),
    ('SLF3.29','Z_6_0',2,1),
    ('SLF3.30','Z_6_2',2,1),
    ('SLF3.31','Z_6_4',2,1),
    ('SLF3.32','Z_6_6',2,1),
    ('SLF3.33','Z_7_1',2,1),
    ('SLF3.34','Z_7_3',2,1),
    ('SLF3.35','Z_7_5',2,1),
    ('SLF3.36','Z_7_7',2,1),
    ('SLF3.37','Z_8_0',2,1),
    ('SLF3.38','Z_8_2',2,1),
    ('SLF3.39','Z_8_4',2,1),
    ('SLF3.40','Z_8_6',2,1),
    ('SLF3.41','Z_8_8',2,1),
    ('SLF3.42','Z_9_1',2,1),
    ('SLF3.43','Z_9_3',2,1),
    ('SLF3.44','Z_9_5',2,1),
    ('SLF3.45','Z_9_7',2,1),
    ('SLF3.46','Z_9_9',2,1),
    ('SLF3.47','Z_10_0',2,1),
    ('SLF3.48','Z_10_2',2,1),
    ('SLF3.49','Z_10_4',2,1),
    ('SLF3.50','Z_10_6',2,1),
    ('SLF3.51','Z_10_8',2,1),
    ('SLF3.52','Z_10_10',2,1),
    ('SLF3.53','Z_11_1',2,1),
    ('SLF3.54','Z_11_3',2,1),
    ('SLF3.55','Z_11_5',2,1),
    ('SLF3.56','Z_11_7',2,1),
    ('SLF3.57','Z_11_9',2,1),
    ('SLF3.58','Z_11_11',2,1),
    ('SLF3.59','Z_12_0',2,1),
    ('SLF3.60','Z_12_2',2,1),
    ('SLF3.61','Z_12_4',2,1),
    ('SLF3.62','Z_12_6',2,1),
    ('SLF3.63','Z_12_8',2,1),
    ('SLF3.64','Z_12_10',2,1),
    ('SLF3.65','Z_12_12',2,1),
    ('SLF3.66','haralick:angular_second_moment',2,1),
    ('SLF3.67','haralick:contrast',2,1),
    ('SLF3.68','haralick:correlation',2,1),
    ('SLF3.69','haralick:sum_of_squares',2,1),
    ('SLF3.70','haralick:inverse_diff_moment',2,1),
    ('SLF3.71','haralick:sum_avg',2,1),
    ('SLF3.72','haralick:sum_var',2,1),
    ('SLF3.73','haralick:sum_entropy',2,1),
    ('SLF3.74','haralick:entropy',2,1),
    ('SLF3.75','haralick:diff_var',2,1),
    ('SLF3.76','haralick:diff_entropy',2,1),
    ('SLF3.77','haralick:info_measure_corr_1',2,1),
    ('SLF3.78','haralick:info_measure_corr_2',2,1),
    ('SLF1.9','edges:area_fraction',2,1),
    ('SLF1.10','edges:homogeneity',2,1),
    ('SLF1.11','edges:direction_maxmin_ratio',2,1),
    ('SLF1.12','edges:direction_maxnextmax_ratio',2,1),
    ('SLF1.13','edges:direction_difference',2,1),
    ('SLF9.1','The number of fluorescent objects in the image',3,1),
    ('SLF9.2','The Euler number of the image',3,1),
    ('SLF9.3','The average object volume',3,1),
    ('SLF9.4','The standard deviation of object volumes',3,1),
    ('SLF9.5','The ratio of the max object volume to min object volume',3,1),
    ('SLF9.6','The average object distance to the protein center of fluorescence (COF)',3,1),
    ('SLF9.7','The standard deviation of object distances from the protein COF',3,1),
    ('SLF9.8','The ratio of the largest to the smallest object to protein COF distance',3,1),
    ('SLF9.9','The average object distance to the COF of the DNA image',3,2),
    ('SLF9.10','The standard deviation of object distances from the COF of the DNA image',3,2),
    ('SLF9.11','The ratio of the largest to the smallest object to DNA COF distance',3,2),
    ('SLF9.12','The distance between the protein COF and the DNA COF',3,2),
    ('SLF9.13','The ratio of the volume occupied by protein to that occupied by DNA',3,2),
    ('SLF9.14','The fraction of the protein fluorescence that co-localizes with DNA',3,2),
    ('SLF9.15','The average horizontal distance of objects to the protein COF',3,1),
    ('SLF9.16','The standard deviation of object horizontal distances from the protein COF',3,1),
    ('SLF9.17','The ratio of the largest to the smallest object to protein COF horizontal distance',3,1),
    ('SLF9.18','The average vertical distance of objects to the protein COF',3,1),
    ('SLF9.19','The standard deviation of object vertical distances from the protein COF',3,1),
    ('SLF9.20','The ratio of the largest to the smallest object to protein COF vertical distance',3,1),
    ('SLF9.21','The average object horizontal distance from the DNA COF',3,2),
    ('SLF9.22','The standard deviation of object horizontal distances from the DNA COF',3,2),
    ('SLF9.23','The ratio of the largest to the smallest object to DNA COF horizontal distance',3,2),
    ('SLF9.24','The average object vertical distance from the DNA COF',3,2),
    ('SLF9.25','The standard deviation of object vertical distances from the DNA COF',3,2),
    ('SLF9.26','The ratio of the largest to the smallest object to DNA COF vertical distance',3,2),
    ('SLF9.27','The horizontal distance between the protein COF and the DNA COF',3,2),
    ('SLF9.28','The signed vertical distance between the protein COF and the DNA COF',3,2),
    ('SLF11.15','The fraction of above threshold pixels that are along an edge',3,1),
    ('SLF11.16','The fraction of fluorescence in above threshold pixels that are along an edge',3,1),
    ('SLF11.17','Average of angular second moment',3,1),
    ('SLF11.18','Average of contrast',3,1),
    ('SLF11.19','Average of correlation',3,1),
    ('SLF11.20','Average of sum of squares of variance',3,1),
    ('SLF11.21','Average of inverse difference moment',3,1),
    ('SLF11.22','Average of sum average',3,1),
    ('SLF11.23','Average of sum variance',3,1),
    ('SLF11.24','Average of sum entropy',3,1),
    ('SLF11.25','Average of entropy',3,1),
    ('SLF11.26','Average of difference variance',3,1),
    ('SLF11.27','Average of difference entropy',3,1),
    ('SLF11.28','Average of info measure of correlation 1',3,1),
    ('SLF11.29','Average of info measure of correlation 2',3,1),
    ('SLF11.30','Range of angular second moment',3,1),
    ('SLF11.31','Range of contrast',3,1),
    ('SLF11.32','Range of correlation',3,1),
    ('SLF11.33','Range of sum of squares of variance',3,1),
    ('SLF11.34','Range of inverse difference moment',3,1),
    ('SLF11.35','Range of sum average',3,1),
    ('SLF11.36','Range of sum variance',3,1),
    ('SLF11.37','Range of sum entropy',3,1),
    ('SLF11.38','Range of entropy',3,1),
    ('SLF11.39','Range of difference variance',3,1),
    ('SLF11.40','Range of difference entropy',3,1),
    ('SLF11.41','Range of info measure of correlation 1',3,1),
    ('SLF11.42','Range of info measure of correlation 2',3,1),
]

def _addfeats(feats):
    global featinfo
    s = len(featinfo)
    featinfo += feats
    e = len(featinfo)
    return range(s,e)
_pftasidxs = _addfeats(pftasinfo())
_overlapidxs = _addfeats(overlapinfo())
_dharidxs = _addfeats(_harprops(1)) + \
            _addfeats(_harprops(2)) + \
            _addfeats(_harprops(3)) + \
            _addfeats(_harprops(4)) + \
            _addfeats(_harprops(5)) + \
            _addfeats(_harprops(6))

_haridxs = range(72, 72+13)
_imgidxs = range(6,20)
_objfieldidxs    = [6, 7, 8, 9, 10,]
_objdnafieldidxs = [6, 7, 8, 9, 10, 18, 19]
_nofidxs = [5]
_edgidxs = range(85,90)
_sklidxs = range(5)

_slf33 = _haridxs + _dharidxs + _objfieldidxs + _edgidxs + _sklidxs + _nofidxs + _pftasidxs
_slf34 = _haridxs + _dharidxs + _objdnafieldidxs + _edgidxs + _sklidxs + _nofidxs + _pftasidxs + _overlapidxs

featuresets = {
    'SLF10' : [103, 113, 102, 93, 91, 90, 104, 110, 111],
    'SLF11' : [90, 91, 92, 93, 94, 95, 96, 97, 104, 105, 106, 107, 108, 109, 118, 119, 120, 121, 122, 123, 124, 125, 126, 127, 128, 129, 130, 131, 132, 133, 134, 135, 136, 137, 138, 139, 140, 141, 142, 143, 144, 145],
    'SLF12' : [8, 80, 25, 5, 77, 82, 29, 85],
    'SLF13' : [8, 80, 19, 25, 5, 77, 29, 82, 85, 74, 11, 16, 87, 7, 43, 30, 2, 66, 18, 76, 6, 56, 83, 10, 72, 14, 4, 13, 86, 75, 73],
    'SLF14' : [90, 91, 92, 93, 94, 95, 96, 97, 104, 105, 106, 107, 108, 109],
    'SLF21' : [0, 1, 2, 3, 4, 8, 9, 10, 72, 73, 74, 75, 76, 77, 78, 79, 80, 81, 82, 83, 84, 85, 86, 87, 88, 89],
    'SLF7' : [6, 7, 8, 9, 10, 11, 12, 13, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79, 80, 81, 82, 83, 84, 85, 86, 87, 88, 89, 5, 0, 1, 2, 3, 4],
    'SLF7DNA' : [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79, 80, 81, 82, 83, 84, 85, 86, 87, 88, 89],
    'SLF8' : [8, 80, 25, 5, 77, 82, 29, 85, 7, 11, 74, 65, 13, 87, 53, 76, 2, 6, 30, 72, 0, 75, 56, 10, 4, 83, 86, 79, 32, 84, 78, 12],
    'SLF9' : [90, 91, 92, 93, 94, 95, 96, 97, 98, 99, 100, 101, 102, 103, 104, 105, 106, 107, 108, 109, 110, 111, 112, 113, 114, 115, 116, 117],
    'SLF33' : _slf33,
    'SLF34' : _slf34,
    '+' : range(len(featinfo)),
}

def _featinfo_for(featset):
    return [featinfo[i] for i in featuresets[featset]]

def get_slf_names(featset):
    return [f[0] for f in _featinfo_for(featset)]
def get_names(featset):
    return [f[1] for f in _featinfo_for(featset)]
def get_channels(featset):
    return [f[3] for f in _featinfo_for(featset)]
