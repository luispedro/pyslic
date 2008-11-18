import numpy
import pyslic
import pyslic.features.tas
from os.path import dirname
basedir=dirname(__file__)

def test_imgfeatcalc():
    Img=pyslic.Image()
    Img.channels[Img.protein_channel]=basedir+'/data/protimg.jp2'
    Img.channels[Img.dna_channel]=basedir+'/data/dnaimg.jp2'
    F=pyslic.computefeatures(Img,'SLF7dna')
    assert F.size == 90

def test_zero_image_wdna():
    img=pyslic.Image()
    img.channeldata[img.protein_channel] = numpy.zeros((100,100))
    img.channeldata[img.dna_channel] = numpy.zeros((100,100))
    img.loaded=True
    img.channels[img.protein_channel]='<special>'
    img.channels[img.dna_channel]='<special>'
    F=pyslic.computefeatures(img,'SLF7dna')
    assert F.size == 90
    assert numpy.isnan(F).sum() == 0

def test_zero_image_nodna():
    img=pyslic.Image()
    img.channeldata[img.protein_channel] = numpy.zeros((100,100))
    img.loaded=True
    img.channels[img.protein_channel]='<special>'
    F=pyslic.computefeatures(img,'mcell')
    assert F.size == 31
    assert numpy.isnan(F).sum() == 0

    F=pyslic.computefeatures(img,['img'])
    assert F.size == 8
    assert numpy.isnan(F).sum() == 0

def test_negative_image():
    img=pyslic.Image()
    img.channeldata[img.protein_channel] = -numpy.ones((10,10))
    img.channeldata[img.dna_channel] = -numpy.ones((10,10))
    img.loaded=True
    img.channels[img.protein_channel]='<special>'
    img.channels[img.dna_channel]='<special>'
    F=pyslic.computefeatures(img,'SLF7dna')
    assert F.size == 90
    assert numpy.isnan(F).sum() == 0

def test_tas():
    img=basedir+'/data/protimg.jp2'
    img=pyslic.image.io.readimg(img)

    IMG=pyslic.Image()
    IMG.channeldata[IMG.protein_channel]=img
    IMG.loaded=True
    IMG.channels[IMG.protein_channel]='<special>'
    assert len(pyslic.features.tas.tas(img)) == 18
    assert len(pyslic.features.tas.pftas(img)) == 18

    assert numpy.all(pyslic.features.tas.tas(img) == pyslic.computefeatures(IMG,['tas']))
