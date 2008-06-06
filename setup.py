from setuptools import setup, find_packages, Extension

convexhull = Extension('imageprocessing/_convexhull', sources = ['imageprocessing/convexhull.cpp'])

magick_libraries= ['Magick++','Wand','Magick'] # This is gotten from magick++-config --libs
readimg = Extension('image.io.readimg', sources = ['image/io/readimg.cpp'], libraries=magick_libraries) 

setup(name='PySLIC',
      version='0.3',
      description='Subcellular Location Image Classifier',
      author='Murphy Lab',
      author_email='murphy@mcu.edu',
      url='http://murphylab.cbi.cmu.edu/',
      packages=find_packages(),
      ext_modules = [convexhull,readimg],
      )


