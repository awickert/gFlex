#! /usr/bin/env python

from setuptools import setup, find_packages

import os

# This is for upload to PyPI
# Should not be necessary on most computers

# Get version from file
import re
VERSIONFILE="gflex/_version.py"
verstrline = open(VERSIONFILE, "rt").read()
VSRE = r"^__version__ = ['\"]([^'\"]*)['\"]"
mo = re.search(VSRE, verstrline, re.M)
if mo:
    __version__ = mo.group(1)
else:
    raise RuntimeError("Unable to find version string in %s." % (VERSIONFILE,))

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
    keywords = ['geophysics', 'geology', 'geodynamics', 'lithosphere', 'isostasy', 'GRASS GIS'],
    classifiers = [],
    url = "https://github.com/awickert/gFlex",
    download_url = "https://github.com/awickert/gFlex/tarball/v"+__version__,
    long_description_content_type='text/markdown',
    long_description = open('README.md', 'r').read(),

    use_2to3=True,
)
