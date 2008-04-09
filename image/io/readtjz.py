import sys
sys.path.append('/home/luispedro/work')
import zipfile
from scipy.misc.pilutil import *
import os
import tempfile
import pyslic
from glob import glob

FILES_PER_TGZ=70

__all__ = ['readtjz_recursive']

def getfileinsidezip(fname,inner):
    '''
    S = getfileinsidezip(zipname,inner)

    Returns the contents of the file inner inside the zip file zipname
    '''
    Z=zipfile.ZipFile(fname)
    return Z.read(inner)

def getimage(fname,inner):
    tmp=tempfile.mkstemp('.jp2')
    tmp2=tempfile.mkstemp('.png')
    S=getfileinsidezip(fname,inner)
    file(tmp[1],'w').write(S)
    os.system("convert %s %s" % (tmp[1],tmp2[1]))
    Img=imread(tmp2[1])
    os.unlink(tmp[1])
    os.unlink(tmp2[1])
    return Img

def readimageinzip(P):
    zip=os.path.dirname(P)
    F=os.path.basename(P)
    return getimage(zip,F)

def parsedir(base):
    Tjzs=glob('%s/*tjz' % base)
    Tjzs.sort()

    images=[]
    for t in Tjzs:
        for i in xrange(FILES_PER_TGZ/2):
            img=pyslic.Image()
            img.set_load_function(readimageinzip)
            p_channel='%s/Stack-%05d' % (t,2*i)
            d_channel='%s/Stack-%05d' % (t,2*i+1)
            img.channels[pyslic.Image.dna_channel]=d_channel
            img.channels[pyslic.Image.protein_channel]=p_channel
            img.label=(t,i)
            images.append(img)
    return images

def readtjz_recursive(base):
    images=[]
    for root,_,_ in os.walk(base):
        images.extend(parsedir(root))
    return images
    

def breaklabel(label):
    T,N=label
    T=os.path.basename(T)
    T=T[:-len('000.flex.tjz')]
    t1,t2=T[:3],T[3:]
    return (int(t1),int(t2),N)

