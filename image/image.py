from scipy.misc.pilutil import *
from numpy import *

__all__ = ['Image', 'setshowimage']

def _open_file_bw(fname):
    """
    Force a file to be B&W
    """
    A=imread(fname)
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

    __slots__ = ['label','scale','regions','channels','channeldata','loaded','__loadfunction']

    def __init__(self):
        '''
        Create an Image.

        Fill the elements in channels to meaningful values.
        '''
        self.label=''
        self.scale = None
        self.channels={}
        self.__loadfunction = _open_file_bw

        self.loaded=False
        self.regions=None
        self.channeldata={}

    def __getstate__(self):
        return (self.label,self.scale,self.channels,self.loaded,self.__loadfunction)

    def __setstate__(self,S):
        self.label,self.scale,self.channels,self.loaded,self.__loadfunction = S

        self.loaded = False
        self.regions = None
        self.channeldata=[]

    def set_load_function(self,f):
        '''
        Sets the file loading function whose signature should be
            string -> numpy array

        The default function does imread() followed by conversion to B&W
            by taking the mean pixel value of the available channels
            (this assumes that the image really is B&W, it was just saved 
            as a colour image)
        '''
        self.__loadfunction = f

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
            if k != self.crop_channel: # Crop is handled like a region
                self.channeldata[k]=self.__loadfunction(v)
        if self.crop_channel in self.channels:
            # These files often need to be fixed:
            self.regions = self.__loadfunction(self.channels[self.crop_channel])
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


    def composite(self,processed = False):
        '''
        composite = img.composite(processed = False)

        Returns a multi colour image containing the different channels
        '''
        self.lazy_load()
        X,Y=self.channeldata[self.protein_channel].shape
        composite=zeros((X,Y,3))
        composite[:,:,1]=self.channeldata[(self.procprotein_channel if processed else self.protein_channel)]
        if self.dna_channel in self.channeldata:
            composite[:,:,0]=self.channeldata[(self.procdna_channel if processed else self.dna_channel)]
        return composite

    def show(self,processed = False):
        '''
        Shows the image composite

        See composite.
        '''
        _showimage(self.composite(processed))

# vim: set ts=4 sts=4 sw=4 expandtab smartindent:
