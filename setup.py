#!/usr/bin/env python3

from setuptools import setup

setup(
    name='mashpit',
    description='A sketch-based surveillance platform',
    scripts=[
        'scripts/create_db.py',
        'scripts/metadata_sra_db.py',
        'scripts/query_against_db.py',
        'scripts/sketch_db.py',
    ],
)
