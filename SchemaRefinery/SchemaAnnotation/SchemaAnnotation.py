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

import argparse

try:
    from SchemaAnnotation.GenbankAnnotations import genbank_annotations, check_genbank_annotations_arguments
    from SchemaAnnotation.ProteomeAnnotations import proteome_annotations, check_proteome_annotations_arguments
    from SchemaAnnotation.AnnotationMerger import annotation_merger
    from SchemaAnnotation.MatchSchemas import check_match_schemas_arguments, match_schemas
    from utils.file_functions import check_and_make_directory

except ModuleNotFoundError:
    from SchemaRefinery.SchemaAnnotation.GenbankAnnotations import genbank_annotations, check_genbank_annotations_arguments
    from SchemaRefinery.SchemaAnnotation.ProteomeAnnotations import proteome_annotations, check_proteome_annotations_arguments
    from SchemaRefinery.SchemaAnnotation.AnnotationMerger import annotation_merger
    from SchemaRefinery.SchemaAnnotation.MatchSchemas import check_match_schemas_arguments, match_schemas
    from SchemaRefinery.utils.file_functions import check_and_make_directory

def check_uniprot_arguments(args_list:list):
    parser = argparse.ArgumentParser(prog='Uniprot Annotations',
                                     description='')
    
    parser.add_argument('-us', '--uniprot_species', type=str, required=True,
                        dest='uniprot_species',
                        help='Uniprot Finder output file for species')

    parser.add_argument('-ug', '--uniprot_genus', type=str, required=False,
                        dest='uniprot_genus',
                        default = '',
                        help='Uniprot Finder output file for genus')

    parser.parse_args(args_list)

def schema_annotation(args):

    error_messages = []
    #check all arguments first
    if 'genbank' in args.annotation_options:

        genbank_arg_list = (["-i", str(args.input_files)] if args.input_files else []) +\
        (["-s" , str(args.schema_directory)] if args.schema_directory else []) +\
        (["-o" , str(args.output_directory)] if args.output_directory else []) +\
        (["-cpu" , str(args.cpu_cores)] if args.cpu_cores else [])

        check_genbank_annotations_arguments(genbank_arg_list)

    if 'proteomes' in args.annotation_options:

        proteomes_arg_list = (["-t", str(args.input_table)] if args.input_table else []) +\
        (["-d", str(args.proteomes_directory)] if args.proteomes_directory else []) +\
        (["-th", str(args.threads)] if args.threads else []) +\
        (["-r" , str(args.retry)] if args.retry else []) +\
        (["-s" , str(args.schema_directory)] if args.schema_directory else []) +\
        (["-o" , str(args.output_directory)] if args.output_directory else []) +\
        (["-cpu" , str(args.cpu_cores)] if args.cpu_cores else [])

        check_proteome_annotations_arguments(proteomes_arg_list)

    if 'uniprot' in args.annotation_options:

        uniprot_arg_list = (["-us", str(args.uniprot_species)] if args.uniprot_species else []) +\
        (["-ug" , str(args.uniprot_genus)] if args.uniprot_genus else [])

        check_uniprot_arguments(uniprot_arg_list)
    
    if 'matchSchemas' in args.annotation_options:

        matchschemas_arg_list = (["-qs", str(args.query_schema)] if args.query_schema else []) +\
        (["-ss", str(args.subject_schema)] if args.subject_schema else []) +\
        (["-oc", str(args.old_schema_columns)] if args.old_schema_columns else []) +\
        (["-ma" , str(args.match_to_add)] if args.match_to_add else []) +\
        (["-o" , str(args.output_directory)] if args.output_directory else []) +\
        (["-bsr" , str(args.blast_score_ratio)] if args.blast_score_ratio else []) +\
        (["-cpu" , str(args.cpu_cores)] if args.cpu_cores else [])

        check_match_schemas_arguments(matchschemas_arg_list)

    check_and_make_directory(args.output_directory)

    # run submodules
    genbank_file = None
    trembl = None
    swiss = None 
    matched_schemas = None
    if 'genbank' in args.annotation_options:
        genbank_file = genbank_annotations(args.input_files, args.schema_directory, args.output_directory, args.cpu_cores)

    if 'proteomes' in args.annotation_options:
        trembl, swiss =proteome_annotations(args.input_table, args.proteomes_directory, args.threads, args.retry, args.schema_directory, args.output_directory, args.cpu_cores)

    if 'matchSchemas' in args.annotation_options:
        matched_schemas = match_schemas(args.query_schema, args.subject_schema, args.output_directory, args.blast_score_ratio, args.cpu_cores)

    annotation_merger(args.uniprot_species, args.uniprot_genus, genbank_file, trembl, swiss, args.match_to_add, args.old_schema_columns, matched_schemas, args.output_directory)