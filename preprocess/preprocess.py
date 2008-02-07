from sys import path
path.append('../')
from image import Image
from basics import fullhistogram
from thresholding import rc
from numpy import *
def preprocess(image,regionid,options):
    def preprocessimg(img):
        img=img.copy()
        cropimg=image.regions
        img = bgsub(img,options)
        img *= (cropimg == regionid)
        T=thresholdfor(img,options)
        img[img < T]=0;
        return img
    image.load()
    img=image.channeldata[Image.protein_channel]
    image.channeldata[Image.procprotein_channel]=preprocessimg(image.channeldata[Image.protein_channel])
    if Image.dna_channel in image.channeldata:
        image.channeldata[Image.procdna_channel]=preprocessimg(image.channeldata[Image.dna_channel])


def thresholdfor(img,options):
    return rc(img)

def bgsub(img,options):
    hist=fullhistogram(img)
    T=argmax(hist[:img.mean()])
    img[img < T]=0
    img=img-T
    return img

# vim: set ts=4 sts=4 sw=4 expandtab smartindent:
