from setuptools import setup, find_packages
import glob
import sys
import os

# VERSION 
VERSION = '0.0.3'

# Get some values from the setup.cfg
try:
    from ConfigParser import ConfigParser
except ImportError:
    from configparser import ConfigParser

conf = ConfigParser()
conf.read(['setup.cfg'])
metadata = dict(conf.items('metadata'))

PACKAGENAME = metadata.get('package_name', 'reduceccd')
DESCRIPTION = metadata.get('description', 'Basic Image CCD reduction using CCDPROC package')
AUTHOR = metadata.get('author', 'Ruben Garcia-Benito')
AUTHOR_EMAIL = metadata.get('author_email', '')
LICENSE = metadata.get('license', 'unknown')
URL = metadata.get('url', 'https://github.com/rgbIAA/reduceccd')

setup(name=PACKAGENAME,
      version=VERSION,
      description=DESCRIPTION,
      install_requires=['astropy', 'numpy', 'scipy', 'ccdproc',
                        'photutils'],
      packages=find_packages(),
      author=AUTHOR,
      author_email=AUTHOR_EMAIL,
      license=LICENSE,
      url=URL
)
