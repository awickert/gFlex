from ez_setup import use_setuptools

use_setuptools ()

from setuptools import setup

setup (name='Flexure',
       version='0.1',
       description='Flexure model',
       author='Andy Wickert',
       author_email='wickert@colorado.edu',
       url='http://csdms.colorado.edu/wiki/Model:Flexure',
       packages=['flexure'],
       scripts=['scripts/flexit.py']
      )
