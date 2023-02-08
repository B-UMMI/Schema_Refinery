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
    from GenbankAnnotations import genbank_annotations
    from ProteomeAnnotations import proteome_annotations
    from AnnotationMerger import annotation_merger

except ModuleNotFoundError:
    from SchemaAnnotation.GenbankAnnotations import genbank_annotations
    from SchemaAnnotation.ProteomeAnnotations import proteome_annotations
    from SchemaAnnotation.AnnotationMerger import annotation_merger, check_uniprot_arguments

def schema_annotation(args):
    if 'genbank' in args.annotation_options:
        genbank_file = genbank_annotations(args.input_files, args.schema_directory, args.output_directory, args.cpu_cores)

    if 'proteomes' in args.annotation_options:
        trembl, swiss =proteome_annotations(args.input_table, args.proteomes_directory, args.threads, args.retry, args.schema_directory, args.output_directory, args.cpu_cores)

    if 'uniprot' in args.annotation_options:
        check_uniprot_arguments(args.species, args.match_to_add, args.matched_schemas)

    annotation_merger(args.species, args.genus, genbank_file, trembl, swiss, args.match_to_add, args.matched_schemas, args.output_directory)