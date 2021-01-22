#!/usr/bin/env python3

from setuptools import setup

setup(
    name='mashpit',
    version='0.7.2',
    url='https://github.com/tongzhouxu/mashpit',
    author='Tongzhou Xu',
    author_email='tongzhou.xu@uga.edu',
    description='A sketch-based surveillance platform',
    packages=['scripts'],
    scripts=[
        'scripts/create_db.py',
        'scripts/metadata_sra_db.py',
        'scripts/query_against_db.py',
        'scripts/sketch_db.py',
    ],
)
