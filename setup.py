#!/usr/bin/env python3

from setuptools import setup

with open("README.md", "r") as rm:
    long_description = rm.read()

setup(
    name='mashpit',
    version='0.8.1',
    url='https://github.com/tongzhouxu/mashpit',
    author='Tongzhou Xu',
    author_email='tongzhou.xu@uga.edu',
    description='A sketch-based surveillance platform',
    long_description=long_description,
    long_description_content_type="text/markdown",
    packages=['mashpit'],
    entry_points={'console_scripts':['mashpit=mashpit.mashpit:main']},
    install_requires=[
        'cython',
        'screed',
        'sourmash',
        'pandas',
        'setuptools',
        'biopython',
        'scipy',
        'python-dotenv'
    ],
)
