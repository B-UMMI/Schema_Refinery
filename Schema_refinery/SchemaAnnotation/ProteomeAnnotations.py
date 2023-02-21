#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Purpose
-------

This script handles everything related with
the creation of TrEMBL and Swiss-Prot annotations.

Code documentation
------------------
"""

import sys

try:
    from SchemaAnnotation.proteome_fetcher import proteome_fetcher
    from SchemaAnnotation.proteome_splitter import proteome_splitter
    from SchemaAnnotation.proteome_matcher import proteome_matcher

except ModuleNotFoundError:
    from Schema_refinery.SchemaAnnotation.proteome_fetcher import proteome_fetcher
    from Schema_refinery.SchemaAnnotation.proteome_splitter import proteome_splitter
    from Schema_refinery.SchemaAnnotation.proteome_matcher import proteome_matcher

def check_proteome_annotations_arguments(input_table, schema_directory):
    necessary_arguments = {
        "input_table": "\tProteome annotations need an input table argument. -t",
        "schema_directory": "\tProteome annotations need a schema directory argument. -s",
    }

    missing_arguments = []

    # Check if all the necessary arguments have been received
    if not input_table:
        missing_arguments.append("input_table")
    
    if not schema_directory:
        missing_arguments.append("schema_directory")

    error_messages = []
    if len(missing_arguments) > 0:
        print("\nError: ")
        for arg in missing_arguments:
            error_messages.append(necessary_arguments[arg])
    
    return error_messages

def proteome_annotations(input_table, proteomes_directory, threads, retry, schema_directory, output_directory, cpu_cores):

    if not proteomes_directory:
        proteome_fetcher(input_table, output_directory, threads, retry)

    tr_file, sp_file, descriptions_file = proteome_splitter(proteomes_directory, output_directory)

    tr_annotations, sp_annotations = proteome_matcher(schema_directory, [tr_file, sp_file, descriptions_file], output_directory, cpu_cores)

    return tr_annotations, sp_annotations
