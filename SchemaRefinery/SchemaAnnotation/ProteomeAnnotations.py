#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Purpose
-------

This sub-module handles everything related with
the creation of TrEMBL and Swiss-Prot annotations.

Code documentation
------------------
"""

import os
import argparse

try:
    from SchemaAnnotation.proteome_fetcher import proteome_fetcher
    from SchemaAnnotation.proteome_splitter import proteome_splitter
    from SchemaAnnotation.proteome_matcher import proteome_matcher
    from utils.file_functions import check_and_make_directory
except ModuleNotFoundError:
    from SchemaRefinery.SchemaAnnotation.proteome_fetcher import proteome_fetcher
    from SchemaRefinery.SchemaAnnotation.proteome_splitter import proteome_splitter
    from SchemaRefinery.SchemaAnnotation.proteome_matcher import proteome_matcher
    from SchemaRefinery.utils.file_functions import check_and_make_directory


def proteome_annotations(input_table, proteomes_directory, threads, retry, schema_directory, output_directory, cpu_cores):

    output_directory = os.path.join(output_directory, "Proteomes")
    check_and_make_directory(output_directory)

    if not proteomes_directory:
        proteomes_directory = proteome_fetcher(input_table, output_directory, threads, retry)

    tr_file, sp_file, descriptions_file = proteome_splitter(proteomes_directory, output_directory)

    tr_annotations, sp_annotations = proteome_matcher(schema_directory, [tr_file, sp_file, descriptions_file], output_directory, cpu_cores)

    return tr_annotations, sp_annotations
