#!/usr/bin/env python

from distutils.core import setup
from setuptools import find_packages

setup(name='py_gaussian',
      version='1.0',
      description='Generating Gaussian16 inputs',
      author='Dominic Willcox',
      author_email='dwillco2@ed.ac.uk',
      packages=find_packages("."),
     )