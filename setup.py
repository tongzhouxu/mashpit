#!/usr/bin/env python3

from setuptools import setup

with open("README.md", "r") as rm:
    long_description = rm.read()

setup(
    name='mashpit',
    version='0.9.6',
    url='https://github.com/tongzhouxu/mashpit',
    author='Tongzhou Xu',
    author_email='tongzhou.xu@uga.edu',
    description='A sketch-based surveillance platform',
    long_description=long_description,
    long_description_content_type="text/markdown",
    packages=['mashpit'],
    include_package_data=True,
    entry_points={'console_scripts':['mashpit=mashpit.mashpit:main']},
    install_requires=[
        'sourmash~=4.6.1',
        'ncbi-datasets-pylib~=14.6.2',
        'pandas',
        'biopython',
        'scikit-bio',
        'tqdm',
        'flask',
        'dask[dataframe]',
    ]
)
