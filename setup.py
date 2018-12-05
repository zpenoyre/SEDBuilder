# to run: python setup.py install

#import numpy as np
try:
    from setuptools import setup
except ImportError:
    from distutils.core import setup
#from distutils.extension import Extension

setup(name='SEDBuilder',
      version='1.0',
      description='Build an SED for an object from VizieR-tabulated photometric data points',
      author='Zephyr Penoyre',
      author_email='zephyrpenoyre@gmail.com',
      url='https://github.com/zpenoyre/SEDBuilder',
      license='MIT',
      packages=['SEDBuilder'],
      #install_requires=['numpy','astroquery','astropy','pathlib','datetime'],
      )
