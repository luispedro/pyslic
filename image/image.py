from scipy.misc.pilutil import *
from numpy import *

__all__ = ['Image', 'setshowimage']

def _open_file(fname):
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

    __slots__ = ['label','scale','regions','channels','channeldata','loaded','__post_load_actions']

    def __init__(self):
        self.label=''
        self.scale = None
        self.regions=None
        self.channels={}
        self.channeldata={}
        self.loaded=False
        self.__post_load_actions = []

    def lazy_load(self):
        '''
        If the image has not been loaded, call load()
        '''
        if not self.loaded:
            self.load()

    def append_post_load(self,action):
        '''
        Register a post load action
        
        action must be callable with a single parameter, the image.
        '''
        self.__post_load_actions.append(action)

    def load(self):
        '''
        Loads the channel data files and regions.

        Calls any post loading actions registered with append_post_load()
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
        for action in self.__post_load_actions:
            action(self)

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
        _showimage(self.composite())

# vim: set ts=4 sts=4 sw=4 expandtab smartindent:
