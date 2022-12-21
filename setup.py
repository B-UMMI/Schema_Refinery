#!/usr/bin/env python

"""The setup script."""

from setuptools import setup

import Schema_refinery

with open('README.rst') as readme_file:
    readme = readme_file.read()

with open('HISTORY.rst') as history_file:
    history = history_file.read()

requirements = [ ]

test_requirements = ['pytest>=3', ]

setup(
    author="UMMI",
    author_email='microbiologia@fm.ul.pt',
    python_requires='>=3.6',
    classifiers=[
        'Development Status :: 2 - Pre-Alpha',
        'Intended Audience :: Developers',
        'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
        'Natural Language :: English',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: 3.8',
    ],
    description="Set of scripts and instructions to refine wg/cgMLST schemas. ",
    install_requires=requirements,
    license="GNU General Public License v3",
    long_description=readme + '\n\n' + history,
    include_package_data=True,
    keywords='schema_refinery',
    name='schema_refinery',
    test_suite='tests',
    tests_require=test_requirements,
    url='https://github.com/B-UMMI/schema_refinery',
    version='0.1.0',
    zip_safe=False,
    
    packages = ['Schema_refinery','Schema_refinery.Download_module'],
    
    entry_points={'console_scripts': ["schema_refinery.py = Schema_refinery.schema_refinery:main",
                                    "SR = Schema_refinery.schema_refinery:main"]}
)
