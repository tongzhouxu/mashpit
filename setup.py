#!/usr/bin/env python3

from setuptools import setup

setup(
    name='mashpit',
    description='A sketch-based surveillance platform',
    packages=['scripts'],
    scripts=[
        'scripts/create_db.py',
        'scripts/metadata_sra_db.py',
        'scripts/query_against_db.py',
        'scripts/sketch_db.py',
    ],
    install_requires=[
        'screed>=1.0.4',
        'sourmash>=3.3.1',
        'setuptools>=47.1.1',
        'biopython>=1.76',
    ],
)
