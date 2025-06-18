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
from setuptools import setup, find_namespace_packages

packages = find_namespace_packages(
  exclude=[]
)

setup(
    packages = packages,
    include_package_data = True
)
