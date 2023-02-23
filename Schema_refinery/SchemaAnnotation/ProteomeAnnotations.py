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
    from Schema_refinery.SchemaAnnotation.proteome_fetcher import proteome_fetcher
    from Schema_refinery.SchemaAnnotation.proteome_splitter import proteome_splitter
    from Schema_refinery.SchemaAnnotation.proteome_matcher import proteome_matcher
    from Schema_refinery.utils.file_functions import check_and_make_directory

def check_proteome_annotations_arguments(args_list:list):
    parser = argparse.ArgumentParser(prog='Proteome Annotations',
                                     description='This sub-module handles everything related with '
                                     'the creation of TrEMBL and Swiss-Prot annotations.')

    parser.add_argument('-s', '--schema-directory', type=str,
                            required=True, dest='schema_directory',
                            help='Path to the schema\'s directory.')

    parser.add_argument('-t', '--input_table', type=str,
                        required=True, dest='input_table',
                        help='TSV file downloaded from UniProt '
                            'that contains list of proteomes.')
    
    parser.add_argument('-d', '--proteomes-directory', type=str, required=False,
                            dest='proteomes_directory',
                            help='Path to the directory with UniProt '
                                'proteomes in Fasta format.')

    parser.add_argument('-th', '--threads', type=int,
                        required=False, default=2,
                        dest='threads',
                        help='Number of threads for concurrent download.')

    parser.add_argument('-r', '--retry', type=int,
                        required=False, dest='retry',
                        default=7,
                        help='Maximum number of retries when a '
                            'download fails.')
    
    parser.add_argument('-cpu', '--cpu-cores', type=int, required=False,
                            dest='cpu_cores',
                            default=1,
                            help='Number of CPU cores to pass to BLAST.')

    parser.add_argument('-o', '--output-directory', type=str,
                        required=True, dest='output_directory',
                        help='Path to the output directory where to save the files.')

    parser.parse_args(args_list)

def proteome_annotations(input_table, proteomes_directory, threads, retry, schema_directory, output_directory, cpu_cores):

    output_directory = os.path.join(output_directory, "Proteomes")
    check_and_make_directory(output_directory)

    if not proteomes_directory:
        proteomes_directory = proteome_fetcher(input_table, output_directory, threads, retry)

    tr_file, sp_file, descriptions_file = proteome_splitter(proteomes_directory, output_directory)

    tr_annotations, sp_annotations = proteome_matcher(schema_directory, [tr_file, sp_file, descriptions_file], output_directory, cpu_cores)

    return tr_annotations, sp_annotations
