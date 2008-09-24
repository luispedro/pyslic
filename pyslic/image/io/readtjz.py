import zipfile
import os
import tempfile
from glob import glob
from ..image import Image
from . import readimg

__all__ = ['readtjz_recursive','readtjz']

def getfileinsidezip(fname,inner):
    '''
    S = getfileinsidezip(zipname,inner)

    Returns the contents of the file inner inside the zip file zipname
    '''
    Z=zipfile.ZipFile(fname)
    S=Z.read(inner)
    Z.close()
    return S

def _getimage(fname,inner):
    '''
    Img = getimage(zipfilename,innername)

    Reads an image inside a zip file
    '''
    S=getfileinsidezip(fname,inner)
    Img=readimg.readimgfromblob(S)
    if len(Img.shape) > 2:
        return Img.mean(2)
    return Img

def readimageinzip(P):
    zip=os.path.dirname(P)
    F=os.path.basename(P)
    return _getimage(zip,F)

def _getlabel(T):
    T=os.path.basename(T)
    if T.endswith('000.flex.tjz'):
        T=T[:-len('000.flex.tjz')]
        t1,t2=T[:3],T[3:]
        return (int(t1),int(t2))
    else:
        return T

def _parsedir(base):
    Tjzs=glob('%s/*000.flex.tjz' % base)
    Tjzs.sort()

    images=[]
    for t in Tjzs:
        images.extend(readtjz(t))
    return images

def readtjz(path):
    '''
    images = readtjz(path)

    Returns all the images inside readtjz.
    '''
    images=[]
    Nstacks=len(filter(lambda inner: inner.startswith('Stack-'),zipfile.ZipFile(path).namelist()))
    for i in xrange(Nstacks/2):
        img=Image()
        img.set_load_function(readimageinzip)
        p_channel='%s/Stack-%05d' % (path,2*i)
        d_channel='%s/Stack-%05d' % (path,2*i+1)
        img.channels[Image.dna_channel]=d_channel
        img.channels[Image.protein_channel]=p_channel
        label=_getlabel(path)
        img.label=label
        img.id=(label,i)
        images.append(img)
    return images

def readtjz_recursive(base):
    '''
    images = readtjz_recursive(basedir)

    Look for all directories below basedir for TJZ files and return the images
    inside them.

    Returns a list of Image objects
    '''
    images=[]
    for root,_,_ in os.walk(base):
        images.extend(_parsedir(root))
    return images
