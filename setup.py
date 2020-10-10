#!/usr/bin/env python3

from setuptools import setup, find_packages

setup(
    name='mashpit',
    description='A sketch-based surveillance platform',
    package_dir={'': 'scripts'},
    packages=find_packages(where='scripts'),
)
