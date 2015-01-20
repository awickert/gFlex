
from ez_setup import use_setuptools
use_setuptools()

from setuptools import setup, find_packages
from setuptools.command.install import install

setup(
    name = "gFlex",
    version = "dev",
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
    description = "Two dimensional plate bending, designed for Earth's lithosphere",
    license = "GPL v3",
    keywords = "geophysics geology geodynamics lithosphere",
    url = "http://csdms.colorado.edu/wiki/Model:GFlex",
    download_url = "https://github.com/awickert/gFlex",
    #long_description=read('README.md'),
)
