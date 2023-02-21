#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Purpose
-------

This module enables the whole process of creating the
Genbank, TrEMBL and Swiss-Prot annotations and merging them
alongside the UniProtFinder and/or match_schemas annotations into a single TSV file.

Code documentation
------------------
"""

import sys

try:
    from SchemaAnnotation.GenbankAnnotations import genbank_annotations, check_genbank_annotations_arguments
    from SchemaAnnotation.ProteomeAnnotations import proteome_annotations, check_proteome_annotations_arguments
    from SchemaAnnotation.AnnotationMerger import annotation_merger
    from SchemaAnnotation.MatchSchemas import check_match_schemas_arguments, match_schemas

except ModuleNotFoundError:
    from Schema_refinery.SchemaAnnotation.GenbankAnnotations import genbank_annotations, check_genbank_annotations_arguments
    from Schema_refinery.SchemaAnnotation.ProteomeAnnotations import proteome_annotations, check_proteome_annotations_arguments
    from Schema_refinery.SchemaAnnotation.AnnotationMerger import annotation_merger
    from Schema_refinery.SchemaAnnotation.MatchSchemas import check_match_schemas_arguments, match_schemas

def check_uniprot_arguments(uniprot_species):

    necessary_arguments = {
        "uniprot_species": "\tUniprot annotations need a species file argument, outputted by Uniprot Finder. -ps",
    }

    missing_arguments = []

    if not uniprot_species:
        missing_arguments.append("uniprot_species")

    error_messages = []
    if len(missing_arguments) > 0:
        print("\nError: ")
        for arg in missing_arguments:
            error_messages.append(necessary_arguments[arg])

    return error_messages

def schema_annotation(args):

    error_messages = []
    #check all arguments first
    if 'genbank' in args.annotation_options:
        error_messages += check_genbank_annotations_arguments(args.input_files, args.schema_directory)

    if 'proteomes' in args.annotation_options:
        error_messages += check_proteome_annotations_arguments(args.input_table, args.schema_directory)

    if 'uniprot' in args.annotation_options:
        error_messages += check_uniprot_arguments(args.uniprot_species)
    
    if 'matchSchemas' in args.annotation_options:
        error_messages += check_match_schemas_arguments(args.query_schema, args.subject_schema, args.old_schema_columns, args.match_to_add)

    if error_messages:
        for message in error_messages:
            print(message)
        sys.exit(0)

    # run submodules
    if 'genbank' in args.annotation_options:
        genbank_file = genbank_annotations(args.input_files, args.schema_directory, args.output_directory, args.cpu_cores)

    if 'proteomes' in args.annotation_options:
        trembl, swiss =proteome_annotations(args.input_table, args.proteomes_directory, args.threads, args.retry, args.schema_directory, args.output_directory, args.cpu_cores)

    if 'matchSchemas' in args.annotation_options:
        matched_schemas = match_schemas(args.query_schema, args.subject_schema, args.output_directory, args.blast_score_ratio, args.cpu_cores)

    annotation_merger(args.uniprot_species, args.uniprot_genus, genbank_file, trembl, swiss, args.match_to_add, args.old_schema_columns, matched_schemas, args.output_directory)