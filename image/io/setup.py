from numpy.distutils.core import setup, Extension

magick_libraries= ['Magick++','Wand','Magick']
readimg = Extension('readimg', sources = ['readimg.cpp'], libraries=magick_libraries) 

setup (name = 'readimg',
       version = '1.0',
       description = 'Read image using ImageMagick',
       ext_modules = [readimg]
       )
