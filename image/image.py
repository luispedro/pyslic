# -*- coding: utf-8 -*-
# Copyright (C) 2008  Murphy Lab
# Carnegie Mellon University
# 
# Written by Luís Pedro Coelho <lpc@cmu.edu>
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
from scipy.misc.pilutil import imshow
from scipy.ndimage import label
import io.readimg
from contextlib import contextmanager

__all__ = ['Image', 'setshowimage','loadedimage']

def _open_file_bw(fname):
    """
    Force a file to be B&W
    """
    A=io.readimg(fname)
    if A.ndim == 3:
        A=A.mean(2)
    return A

_showimage=imshow
def setshowimage(f):
    '''
    Set the function that shows an image.

    By default, uses scipy.misc.pilutil.imshow.
    '''
    global _showimage
    _showimage = f


@contextmanager
def loadedimage(img):
    '''
    with loadedimage(img) as img:
        ...

    On entry, loads the image if needed (equivalent to img.lazy_load()).
    On exit, restores the previous state.
    '''
    wasloaded=img.loaded
    img.lazy_load()
    yield img
    if not wasloaded:
        img.unload()

class Image(object):
    """
    Represents a multi-channel image.

    Images are loaded in a lazy fashion.
    """
    dna_channel='dna'
    protein_channel='protein'
    autofluorescence_channel='autofluorescence'
    crop_channel='crop'

    procdna_channel='procdna'
    procprotein_channel='procprotein'
    residualprotein_channel='resprotein'

    __slots__ = [
                'attributes',
                'label',
                'id',
                'scale',
                'regions',
                'channels',
                'channeldata',
                'loaded',
                'load_function',
                ]

    def __init__(self):
        '''
        Create an Image.

        Fill the elements in channels to meaningful values.
        '''
        self.label=''
        self.id=None
        self.scale = None
        self.channels={}
        self.load_function = _open_file_bw

        self.loaded=False
        self.regions=None
        self.channeldata={}
        self.attributes={}

    def __getstate__(self):
        return (self.attributes,self.label,self.id,self.scale,self.channels,self.loaded,self.load_function)

    def __setstate__(self,S):
        if len(S) == 7:
            self.attributes,self.label,self.id,self.scale,self.channels,self.loaded,self.load_function = S
        else:
            self.label,self.id,self.scale,self.channels,self.loaded,self.load_function = S
            self.attributes = {}

        self.loaded = False
        self.regions = None
        self.channeldata={}

    def set_load_function(self,f):
        '''
        Sets the file loading function whose signature should be
            string -> numpy array

        The default function does imread() followed by conversion to B&W
            by taking the mean pixel value of the available channels
            (this assumes that the image really is B&W, it was just saved 
            as a colour image)
        '''
        self.load_function = f

    def lazy_load(self):
        '''
        If the image has not been loaded, call load()
        '''
        if not self.loaded:
            self.load()

    def load(self):
        '''
        Loads the channel data files and regions.

        Calls any post loading actions registered with append_post_load()
        '''
        for k,v in self.channels.items():
            if k != self.crop_channel: # Crop is called "regions"
                self.channeldata[k]=self.load_function(v)
        if self.crop_channel in self.channels:
            self.regions = self.load_function(self.channels[self.crop_channel])
            # These files often need to be fixed 
            self.regions,_ = label(self.regions)
        self.loaded = True

    def unload(self):
        '''
        Unloads the channel data
        '''
        self.channeldata={}
        self.regions=None
        self.loaded=False


    def composite(self,processed = False):
        '''
        composite = img.composite(processed = False)

        Returns a multi colour image containing the different channels
        '''
        self.lazy_load()
        X,Y=self.channeldata[self.protein_channel].shape
        composite=numpy.zeros((X,Y,3))
        composite[:,:,1]=self.channeldata[(self.procprotein_channel if processed else self.protein_channel)]
        if self.dna_channel in self.channeldata:
            composite[:,:,0]=self.channeldata[(self.procdna_channel if processed else self.dna_channel)]
        if self.autofluorescence_channel in self.channeldata:
            composite[:,:,2]=self.channeldata[self.autofluorescence_channel]
        return composite

    def show(self,processed = False):
        '''
        Shows the image composite

        See composite.
        '''
        _showimage(self.composite(processed))

# vim: set ts=4 sts=4 sw=4 expandtab smartindent:
