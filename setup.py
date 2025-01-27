#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Purpose
-------
This script installs the schema_refinery package in the local files, giving the
option to call schema_rifinery in the command line by typing SR or schema_refinery.

Code documentation
------------------
"""
from SchemaRefinery.utils import constants as ct
from setuptools import setup, find_namespace_packages

packages = find_namespace_packages(
  exclude=[]
)

with open('README.rst',encoding='utf-8') as readme_file:
    readme = readme_file.read()

with open('HISTORY.rst',encoding='utf-8') as history_file:
    history = history_file.read()

test_requirements = ['pytest>=3', ]

setup(
    author="UMMI",
    author_email='imm-bioinfo@medicina.ulisboa.pt',
    python_requires='>=3.9',
    install_requires = ct.DEPENDENCIES_VERSION,
    classifiers=[
        'Development Status :: 3 - Alpha',
        'Intended Audience :: Developers',
        'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
        'Natural Language :: English',
        'Operating System :: Unix',
        'Environment :: Console',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.9',
        'Programming Language :: Python :: 3.10',
        'Programming Language :: Python :: 3.11',
        'Programming Language :: Python :: 3.12',
        'Programming Language :: Python :: 3.13',
        'Programming Language :: Python :: Implementation :: CPython',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
    ],
    description="Tool to refine cg/wgMLST Schemas",
    license="GNU General Public License v3",
    long_description=readme + '\n\n' + history,
    include_package_data=True,
    keywords='schema_refinery',
    name='schema_refinery',
    test_suite='tests',
    tests_require=test_requirements,
    url='https://github.com/B-UMMI/schema_refinery',
    version=ct.VERSION,
    zip_safe=False,

    packages = packages,

    entry_points={'console_scripts': ["SchemaRefinery = SchemaRefinery.schema_refinery:entry_point",
                                    "SR = SchemaRefinery.schema_refinery:entry_point"]}
)
