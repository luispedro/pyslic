from distutils.core import setup, Extension

convexhull = Extension('_convexhull', sources = ['convexhull.cpp'])

setup (name = 'Convex Hull',
       version = '1.0',
       description = 'Support package for feature calculation',
       ext_modules = [convexhull])
