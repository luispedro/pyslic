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
import numpy as np
from numpy import *
from scipy import ndimage
from mahotas.edge import sobel
import math

def edgefeatures(protproc):
    """
    values = edgefeatures(protproc)
       where protproc contains the pre-processed fluorescence image.
       Pre-processed means that the image has been cropped and had
      Features calculated include:
      (The following feature descriptions were added here by T. Zhao
      according to the reference
      "R. F. Murphy, M. Velliste, and G. Porreca (2003) Robust Numerical
      Features for Description and Classification of Subcellular Location
      Patterns in Fluorescence Microscope Images. J. VLSI Sig. Proc. 35:
      311-321.")
      1. Fraction of above-threshold pixels along edge
      2. Measure of edge gradient intensity homogenerity
      3. Measure of edge direction homogenerity 1
      4. Measure of edge direction homogenerity 2
      5. Measure of edge direction difference

    07 Mar 99 - M.V. Boland

    Corrections by M. Velliste 1/20/02
    1) corrected the way non-zero edge pixels are selected (see below).
    2) corrected the homogeneity feature so that it is based on a
       histogram of edge magnitudes (see below).
    3) Made some comments below about differences between Matlab 5 &
       6, and provided a fix for version 5 so that it would give
       consistent results.

    M.Velliste June 2, 2002: added SLF names
    Ported to Python by Luis Pedro Coelho
    """

    binimg = (protproc > 0)
    edges = sobel(protproc)
    A = edges.sum()/binimg.sum()
    #A = bwarea(edge(imageproc,'canny',[]))/bwarea(im2bw(imageproc)) ;

    # Directional edge filters
    N = np.array([
            [ 1, 1, 1],
            [ 0, 0, 0],
            [-1,-1,-1]],
            np.float_)
    W = np.array([
            [ 1, 0,-1],
            [ 1, 0,-1],
            [ 1, 0,-1]],
            np.float_)

    # Calculation of the gradient from two orthogonal directions
    protproc = asarray(protproc.copy(), np.float_)
    iprocN = ndimage.convolve(protproc,N)
    iprocW = ndimage.convolve(protproc,W)

    # Calculate the magnitude and direction of the gradient
    iprocmag = np.sqrt(iprocN**2. + iprocW**2.)
    iproctheta = np.arctan2(iprocN, iprocW)

    # Change by MV:
    # Identify pixels in iprocmag that are not 0
    #  (i.e. iprocN and iprocW were both 0). Before
    # this was incorrectly based on identifying non-zero
    # pixels in iproctheta, which does remove zero-magnitude
    # edges, but it also removes edges that face exactly east.
    v = iproctheta.ravel()
    v_mag = iprocmag.ravel()
    v = v[v_mag > 0]
    v_mag = v_mag[v_mag > 0]

    if v.size == 0:
        return zeros(5)

    # Histogram the gradient directions
    h,_ = histogram(v,8)

    # max/min ratio
    maxidx=argmax(h)
    hmax=h[maxidx]
    hmin=h.min()

    if hmin > 0:
        maxminratio = hmax/hmin
    else:
        maxminratio = 0


    # Difference between bins of histogram at angle and angle+pi
    #  In general, objects have an equal number of pixels at an angle
    #  and that angle+pi. The differences are normalized to the sum of
    #  the two directions.
    diff = abs(h[:4]-h[4:])/abs(h[:4]+h[4:])
    diff[abs(h[:4]-h[4:])==0] = 0


    h[maxidx] = 0
    hnextmax = h.max()
    maxnextmaxratio=hmax/hnextmax

    sumdiff=diff.sum()

    # Measure of edge homogeneity - what fraction of edge pixels are in
    #  the first two bins of the histogram.
    # Change by MV: Made it be based on edge magnitude histogram. Was
    # incorrectly based on edge direction histogram before.
    h_mag,_ = histogram(v_mag,4)
    homogeneity = h_mag[0]/h_mag.sum()

    return array([A,homogeneity,maxminratio,maxnextmaxratio,sumdiff])

edgefeatures.names = ['edges:area_fraction', 'edges:homogeneity', 'edges:direction_maxmin_ratio', 'edges:direction_maxnextmax_ratio', 'edges:direction_difference']
edgefeatures.slfnames = ['SLF1.9', 'SLF1.10', 'SLF1.11', 'SLF1.12', 'SLF1.13']
# vim: set ts=4 sts=4 sw=4 expandtab smartindent:
