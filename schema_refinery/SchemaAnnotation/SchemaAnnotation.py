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

try:
    from SchemaAnnotation.GenbankAnnotations import genbank_annotations
    from SchemaAnnotation.ProteomeAnnotations import proteome_annotations

except ModuleNotFoundError:
    from SchemaAnnotation.GenbankAnnotations import genbank_annotations
    from SchemaAnnotation.ProteomeAnnotations import proteome_annotations

def schema_annotation(args):
    if 'genbank' in args.annotation_options:
        genbank_annotations(args.input_files, args.schema_directory, args.output_directory, args.cpu_cores)

    if 'proteomes' in args.annotation_options:
        proteome_annotations(args.input_table, args.proteomes_directory, args.threads, args.retry, args.schema_directory, args.output_directory, args.cpu_cores)

   # TODO run annotation_merger