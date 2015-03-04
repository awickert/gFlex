#! /usr/bin/env python

from ez_setup import use_setuptools
use_setuptools()

from setuptools import setup, find_packages
from setuptools.command.install import install

import os

# This is for upload to PyPI
# Should not be necessary on most computers
try:
    long_description = open('README.txt').read()
except:
    long_description = "see README.md"
setup(
    name = "gFlex",
    version = "0.8.1",
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
    download_url = "https://github.com/awickert/gFlex/tarball/v0.8.1",
    long_description = long_description,
)
