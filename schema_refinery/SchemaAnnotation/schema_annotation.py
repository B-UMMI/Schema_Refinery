#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Purpose
-------

This module enables the whole process of creating the
Genbank, TrEMBL and Swiss-Prot annotations and merging them
alongside the UniProtFinder annotations into a single TSV file.

Code documentation
------------------
"""

import os

from genbank_annotations import genbankAnnotations
from proteome_annotations import proteomeAnnotations

# GENBANK ANNOTATIONS
# parser.add_argument('-i', type=str, required=True,
#                         dest='input_files',
#                         help='Path to the directory that contains '
#                              'Genbank files with annotations to '
#                              'extract.')

# PROTEOME ANNOTATIONS
# parser.add_argument('-t', '--input_table', type=str,
#                     required=True, dest='input_table',
#                     help='TSV file downloaded from UniProt '
#                          'that contains list of proteomes.')

# parser.add_argument('-i', type=str, required=True,
#                         dest='proteomes_directory',
#                         help='Path to directory with UniProt '
#                              'proteomes in Fasta format.')

# parser.add_argument('-th', '--threads', type=int,
#                     required=False, default=2,
#                     dest='threads',
#                     help='Number of threads for concurrent download.')

# parser.add_argument('-r', '--retry', type=int,
#                     required=False, dest='retry',
#                     default=7,
#                     help='Maximum number of retries when a '
#                          'download fails.')

# REPEATED

# parser.add_argument('-s', '--schema-directory', type=str,
#                         required=True, dest='schema_directory',
#                         help='Path to the schema\'s directory.')

# parser.add_argument('-o', '--output-directory', type=str,
#                         required=True, dest='output_directory',
#                         help='Path to the output directory.')

# parser.add_argument('-cpu', '--cpu-cores', type=int, required=False,
#                         dest='cpu_cores',
#                         default=1,
#                         help='Number of CPU cores to pass to BLAST.')

def main(input_files:str, input_table:str, proteomes_directory:str, schema_directory:str, threads:int, cpu_cores:int, retry:int, output_directory:str):
    proteomes_output_directory = os.path.join(output_directory, "Proteome_Annotations_Output")
    genbank_output_directory = os.path.join(output_directory, "Genbank_Annotations_Output")

    try:
        genbankAnnotations(input_files, schema_directory, genbank_output_directory, cpu_cores)
        proteomeAnnotations(input_table, proteomes_directory, schema_directory, threads, retry, cpu_cores, proteomes_output_directory)

        # TODO run annotation_merger


    except:
        print("Something went wrong")
