# -*- coding: utf-8 -*-
# Copyright (C) 2008  Murphy Lab
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

from __future__ import division
import numpy
from scipy.ndimage import label
from contextlib import contextmanager
from mahotas.freeimage import imread

__all__ = ['Image', 'setshowimage','loadedimage']

def _open_file_bw(fname):
    """
    Force a file to be B&W
    """
    A = imread(fname)
    if A.ndim == 3:
        A = A.max(2)
    return A

def _showimage(img):
    from pylab import imshow
    #from scipy.misc.pilutil import imshow
    imshow(img)

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

    Class variables:
        * label: a label
        * id: an id.
            'label' & 'id' are for the user (i.e., class Image will not use them/change them). Generally, 
            in a given collection 'id' is unique, but 'label' is shared across many images in the same class.

        * channels: this is a dictionary from channel-id to filename

        * channeldata: Once the image is loaded, this is a dictionary from channel-id to image array

        * loaded: Boolean, whether the image has been loaded.
        * load_function: a function that takes a filename and returns an array.
        * post_load: a list of functions that are called post-loading. They are called with self as their single argument.
            This is useful for implementing normalisation, for example.
                
        * scale: scale of the image in microns/pixel.
        * regions: saves the segmentation of the image as a labeled map of regions.


    Pickling
    --------

    Images can be pickled, but image data *does not* go with the pickling. Pickling an image is always the same as pickling its
    unloaded version.
    """

    procdna_channel='procdna'
    procprotein_channel='procprotein'
    residualprotein_channel='resprotein'


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
        self.post_load=[]
        self.temp = {}

    def __getstate__(self):
        copy = self.__dict__.copy()
        del copy['channeldata']
        del copy['loaded']
        del copy['regions']
        del copy['temp']
        items = copy.items()
        items.sort(key=lambda it: it[0])
        return items

    def __setstate__(self, state):
        self.loaded = False
        self.channeldata = {}
        self.regions = None
        self.temp = {}
        for k,v in state:
            self.__dict__[k] = v

    def __repr__(self):
        '''Implement repr() operator'''
        return 'Image( %s )' % repr(self.channels)

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
            if k != 'crop': # Crop is called "regions"
                if type(v) == list:
                    self.channeldata[k]=[]
                    for f in v:
                        self.channeldata[k].append(self.load_function(f))
                    self.channeldata[k]=numpy.array(self.channeldata[k])
                else:
                    self.channeldata[k]=self.load_function(v)
        if 'crop' in self.channels:
            self.regions = self.load_function(self.channels['crop'])
            # These files often need to be fixed 
            self.regions,_ = label(self.regions)
        self.loaded = True
        for post in self.post_load:
            post(self)

    def unload(self):
        '''
        Unloads the channel data
        '''
        self.channeldata={}
        self.regions=None
        self.loaded=False

    def nr_slices(self):
        '''
        nr_slices = img.nr_slices()

        Returns the number of z slices.
        '''
        if type(self.channels['protein']) != list:
            return 1
        return len(self.channels['protein'])

    def get(self,channelid):
        '''
        ch = img.get(channelid)

        Equivalent of
            img.lazy_load()
            ch = img.channeldata[channelid]
        '''
        self.lazy_load()
        return self.channeldata[channelid]

    def composite(self, idx = 0, processed = False):
        '''
        composite = img.composite(idx = 0, processed = False)

        Returns a multi colour image containing the different channels
        '''
        self.lazy_load()
        def getchannel(channel):
            if type(self.channels['protein']) == list:
                orig=self.channeldata[channel][idx,:,:]
            else:
                orig=self.channeldata[channel]
            if orig.ptp():
                return numpy.array( (orig.astype(numpy.float)-orig.min()) * 255./orig.ptp(), numpy.uint8 )
            return orig
        prot = getchannel('protein')
        X,Y=prot.shape
        composite=numpy.zeros((X,Y,3),numpy.uint8)
        composite[:,:,1]=getchannel('procprotein' if processed else 'protein')
        if 'dna' in self.channeldata:
            composite[:,:,0] = getchannel('procdna' if processed else 'dna')
        if 'autofluorescence' in self.channeldata:
            composite[:,:,2] = getchannel('autofluorescence')
        return composite

    def show(self, idx = 0, processed = False):
        '''
        Shows the image composite

        See composite.
        '''
        _showimage(self.composite(idx,processed))

# vim: set ts=4 sts=4 sw=4 expandtab smartindent:
