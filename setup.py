#!/usr/bin/env python3.4

from setuptools import setup,find_packages

setup(name='BioTools',
      version='1.0',
      author='Jacob Brady',
      author_email='jacob.brady0449@gmail.com',
      url='https://github.com/j-brady/BioTools',

      packages=find_packages(),
      description='Python scripts for biochemists and NMR spectroscopists',
      
      install_requires = ['numpy','matplotlib','BioPython']
)
