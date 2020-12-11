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
        'cPython>=0.0.6',
        'screed>=1.0.4',
        'sourmash>=3.3.1',
        'pandas>=1.0.3',
        'setuptools>=47.1.1',
        'biopython>=1.76',
        'numpy>=1.16.5',
    ],
)
