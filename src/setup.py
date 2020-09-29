#!/usr/bin/python

from distutils.core import setup
from distutils.core import Extension
import os


module=Extension('eagle', 
                 include_dirs = [os.environ['PYTHONPATH']+'/../python3.5/site-packages/numpy/core/include'],
                 libraries = ['hdf5'],
                 sources = ['readEagle.cpp']
                 )

setup(name = 'Eagle',
      version = '1.0',
      description='Read EAGLE files in serial or parallel mode.',
      ext_modules = [module],
      author='Matthieu Schaller',
      author_email='matthieu.schaller@durham.ac.uk',
      url='http://eagle.strw.leidenuniv.nl/wiki/doku.php?id=eagle:documentation:reading_python',
      long_description=''' See website. '''
      )
