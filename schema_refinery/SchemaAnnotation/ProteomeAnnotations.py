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
    from proteome_fetcher import proteome_fetcher
    from proteome_splitter import proteome_splitter
    from proteome_matcher import proteome_matcher

except ModuleNotFoundError:
    from SchemaAnnotation.proteome_fetcher import proteome_fetcher
    from SchemaAnnotation.proteome_splitter import proteome_splitter
    from SchemaAnnotation.proteome_matcher import proteome_matcher

def proteome_annotations(input_table, proteomes_directory, threads, retry, schema_directory, output_directory, cpu_cores):

    necessary_arguments = {
        "input_table": "\tProteome annotations need an input table argument. -t",
        "proteomes_directory": "\tProteome annotations need a proteomes directory argument. -d",
        "schema_directory": "\tProteome annotations need a schema directory argument. -s",
        "output_directory": "\tProteome annotations need an output directory argument. -o"
    }

    missing_arguments = []

    # Check if all the necessary arguments have been received
    if not input_table:
        missing_arguments.append("input_table")
    
    if not proteomes_directory:
        missing_arguments.append("proteomes_directory")
    
    if not schema_directory:
        missing_arguments.append("schema_directory")
    
    if not output_directory:
        missing_arguments.append("output_directory")

    if len(missing_arguments > 0):
        print("\nError: ")
        for arg in missing_arguments:
            print(necessary_arguments[arg])
        sys.exit(0)

    # proteomes are downloaded to proteomes_directory
    proteome_fetcher(input_table, proteomes_directory, threads, retry)

    tr_file, sp_file, descriptions_file = proteome_splitter(proteomes_directory, output_directory)

    tr_annotations, sp_annotations = proteome_matcher(schema_directory, [tr_file, sp_file, descriptions_file], output_directory, cpu_cores)

    return tr_annotations, sp_annotations
