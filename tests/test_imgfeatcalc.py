import numpy
import pyslic
import pyslic.features.tas
from os.path import dirname
basedir=dirname(__file__)

def test_imgfeatcalc():
    Img=pyslic.Image()
    Img.channels['protein']=basedir+'/data/protimg.jp2'
    Img.channels['dna']=basedir+'/data/dnaimg.jp2'
    F=pyslic.computefeatures(Img,'SLF7dna')
    assert F.size == 90

def test_zero_image_wdna():
    img=pyslic.Image()
    img.channeldata['protein'] = numpy.zeros((100,100))
    img.channeldata['dna'] = numpy.zeros((100,100))
    img.loaded=True
    img.channels['protein']='<special>'
    img.channels['dna']='<special>'
    F=pyslic.computefeatures(img,'SLF7dna')
    assert F.size == 90
    assert numpy.isnan(F).sum() == 0

def test_zero_image_nodna():
    img=pyslic.Image()
    img.channeldata['protein'] = numpy.zeros((100,100))
    img.loaded=True
    img.channels['protein']='<special>'
    F=pyslic.computefeatures(img,'mcell')
    assert F.size == 31
    assert numpy.isnan(F).sum() == 0

    F=pyslic.computefeatures(img,['img'])
    assert F.size == 8
    assert numpy.isnan(F).sum() == 0

def test_negative_image():
    img=pyslic.Image()
    img.channeldata['protein'] = -numpy.ones((10,10))
    img.channeldata['dna'] = -numpy.ones((10,10))
    img.loaded=True
    img.channels['protein']='<special>'
    img.channels['dna']='<special>'
    F=pyslic.computefeatures(img,'SLF7dna')
    assert F.size == 90
    assert numpy.isnan(F).sum() == 0

def test_tas():
    img=basedir+'/data/protimg.jp2'
    img=pyslic.image.io.readimg(img)

    IMG=pyslic.Image()
    IMG.channeldata['protein']=img
    IMG.loaded=True
    IMG.channels['protein']='<special>'
    assert len(pyslic.features.tas.tas(img)) == 2*3*9
    assert len(pyslic.features.tas.pftas(img)) == 2*3*9

    assert numpy.all(pyslic.features.tas.tas(img) == pyslic.computefeatures(IMG,['tas']))

def test_imgfeaturesfield():
    img = basedir+'/data/protimg.jp2'
    img = pyslic.image.io.readimg(img)
    assert len(pyslic.features.imgfeatures.imgfeaturesdna(img,None,True)) < len(pyslic.features.imgfeatures.imgfeaturesdna(img,None,False))

def test_featinfolen():
    imgdata = basedir+'/data/protimg.jp2'
    imgdata = pyslic.image.io.readimg(imgdata)
    img = pyslic.Image()
    img.channeldata['protein'] = imgdata
    img.loaded = True
    img.channels['protein']='<special>'
    assert len(pyslic.computefeatures(img,'SLF33')) == len(pyslic.features.featinfo.get_names('SLF33'))
