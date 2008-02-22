from scipy.misc.pilutil import *
from numpy import *

__all__ = ['Image']

def _open_file(fname):
    """
    Force a file to be B&W
    """
    A=imread(fname)
    if A.ndim == 3:
        A=A.mean(2)
    return A

class Image(object):
    """
    class Image(object)

    Represents a multi-channel image.
    """
    dna_channel='dna'
    protein_channel='protein'
    autofluorescence_channel='autofluorescence'
    crop_channel='crop'

    procdna_channel='procdna'
    procprotein_channel='procprotein'
    residualprotein_channel='resprotein'

    __slots__ = ['label','scale','regions','channels','channeldata','loaded']
    def __init__(self):
        self.label=''
        self.scale = None
        self.regions=None
        self.channels={}
        self.channeldata={}
        self.loaded=False

    def lazy_load(self):
        '''
        If the image has not been loaded, call load()
        '''
        if not self.loaded:
            self.load()

    def load(self):
        '''
        Loads the channel data files and regions
        '''
        for k,v in self.channels.items():
            if k != self.crop_channel: # Crop is handled like a region
                self.channeldata[k]=_open_file(v)
        if self.crop_channel in self.channels:
            # These files often need to be fixed:
            self.regions = _open_file(self.channels[self.crop_channel])
            if self.regions.max() == 255:
                self.regions[self.regions == 255] = 1
        self.loaded = True

    def unload(self):
        '''
        Unloads the channel data
        '''
        self.channeldata={}
        self.regions=None
        self.loaded=False


    def composite(self):
        '''
        composite = img.composite()

        Returns a multi colour image containing the different channels
        '''
        self.lazy_load()
        X,Y=self.channeldata[self.protein_channel].shape
        composite=zeros((X,Y,3))
        composite[:,:,1]=self.channeldata[self.protein_channel]
        if self.dna_channel in self.channeldata:
            composite[:,:,0]=self.channeldata[self.dna_channel]
        return composite

    def show(self):
        '''
        Shows the image composite

        See composite.
        '''
        imshow(self.composite())

# vim: set ts=4 sts=4 sw=4 expandtab smartindent:
