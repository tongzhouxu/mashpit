#!/usr/bin/env python3

from setuptools import setup

with open("README.md", "r") as rm:
    long_description = rm.read()

setup(
    name='mashpit',
    version='0.9.2',
    url='https://github.com/tongzhouxu/mashpit',
    author='Tongzhou Xu',
    author_email='tongzhou.xu@uga.edu',
    description='A sketch-based surveillance platform',
    long_description=long_description,
    long_description_content_type="text/markdown",
    packages=['mashpit'],
    entry_points={'console_scripts':['mashpit=mashpit.mashpit:main']},
    install_requires=[
        'numpy~=1.19.5',
        'screed~=1.0.5',
        'sourmash~=4.2.2',
        'pandas~=1.1.5',
        'biopython~=1.78',
        'scipy~=1.7.3',
        'python-dotenv',
        'setuptools==65.5.1'
    ]
)
