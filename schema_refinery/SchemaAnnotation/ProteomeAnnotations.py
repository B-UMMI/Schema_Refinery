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

    # Check if all the necessary arguments have been received
    if not input_table:
        sys.exit("\nError: Proteome annotations need an input table argument. -t")
    
    if not proteomes_directory:
        sys.exit("\nError: Proteome annotations need a proteomes directory argument. -d")

    if not schema_directory:
        sys.exit("\nError: Proteome annotations need a schema directory argument. -s")
    
    if not output_directory:
        sys.exit("\nError: Proteome annotations need an output directory argument. -o")

    # proteomes are downloaded to proteomes_directory
    # proteome_fetcher(args.input_table, args.proteomes_directory, args.threads, args.retry)

    # tr_file, sp_file, descriptions_file = proteome_splitter(args.proteomes_directory, args.output_directory)

    # proteomeMatcher(args.schema_directory, [tr_file, sp_file, descriptions_file], args.output_directory, args.cpu_cores)
