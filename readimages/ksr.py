import re
import os.path
import os
import sys
from ..image import Image

__all__ = ['read_ksr_dir']

def read_ksr_dir(dir):
    Files=os.listdir(dir)
    pat=re.compile('KSR_.*t([0-9]+)([A-H][0-9]{1,2})f([0-9]+)d([0-9])\.(tif|TIF)')
    channelcode = { 1 : Image.dna_channel , 2 : Image.protein_channel, 3 : Image.autofluorescence_channel } 
    images={}
    for f in Files:
        m=pat.match(f)
        if not m:
            print "Don't know how to process", f
            continue
        T,Well,Field,Channel,_=m.groups()
        Field=int(Field)
        img=images.get((Well,Field),None)
        if img is None:
            img = Image()
            images[(Well,Field)]=img
        imgs.channels[channelcode[int(Channel)]]=os.path.join(dir,f)
    return images


# vim: set ts=4 sts=4 sw=4 expandtab smartindent:
