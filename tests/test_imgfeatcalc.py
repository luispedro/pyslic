import numpy
import pyslic
from os.path import dirname
basedir=dirname(__file__)

def test_imgfeatcalc():
    Img=pyslic.Image()
    Img.channels[Img.protein_channel]=basedir+'/data/protimg.jp2'
    Img.channels[Img.dna_channel]=basedir+'/data/dnaimg.jp2'
    F=pyslic.computefeatures(Img,'SLF7dna')
    assert F.size == 90

def test_zero_image():
    img=pyslic.Image()
    img.channeldata[img.protein_channel] = numpy.zeros((10,10))
    img.loaded=True
    img.channels[img.protein_channel]='<special>'
    F=pyslic.computefeatures(img,'SLF7dna')
    assert F.size == 90
    assert numpy.isnan(F).sum() == 0

