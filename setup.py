#! /usr/bin/env python

from ez_setup import use_setuptools
use_setuptools()

from setuptools import setup, find_packages
from setuptools.command.install import install

import os

# This is for upload to PyPI
# Should not be necessary on most computers

import re
VERSIONFILE="gflex/_version.py"
verstrline = open(VERSIONFILE, "rt").read()
VSRE = r"^__version__ = ['\"]([^'\"]*)['\"]"
mo = re.search(VSRE, verstrline, re.M)
if mo:
    __version__ = mo.group(1)
else:
    raise RuntimeError("Unable to find version string in %s." % (VERSIONFILE,))

try:
    long_description = open('README.txt').read()
except:
    long_description = "see README.md"
setup(
    name = "gFlex",
    version = __version__,
    packages = find_packages(exclude="tests"),
    entry_points = {
      'console_scripts': ['gflex = gflex:main']
      },

    package_data = { 
      '': ['*.md']
      },

    # metadata for upload to PyPI
    author = "Andrew D. Wickert",
    author_email = "awickert@umn.edu",
    description = "One- and two-dimensional plate bending, designed for Earth's lithosphere",
    license = "GPL v3",
    keywords = ['geophysics', 'geology', 'geodynamics', 'lithosphere', 'isostasy'],
    classifiers = [],
    url = ["https://github.com/awickert/gFlex", "http://csdms.colorado.edu/wiki/Model:gFlex"],
    download_url = "https://github.com/awickert/gFlex/tarball/v"+__version__,
    long_description = long_description,
)
