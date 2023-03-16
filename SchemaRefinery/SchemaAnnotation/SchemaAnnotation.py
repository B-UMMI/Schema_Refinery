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


try:
    from SchemaAnnotation.GenbankAnnotations import genbank_annotations
    from SchemaAnnotation.ProteomeAnnotations import proteome_annotations
    from SchemaAnnotation.AnnotationMerger import annotation_merger
    from SchemaAnnotation.MatchSchemas import match_schemas
    from utils.file_functions import check_and_make_directory
except ModuleNotFoundError:
    from SchemaRefinery.SchemaAnnotation.GenbankAnnotations import genbank_annotations
    from SchemaRefinery.SchemaAnnotation.ProteomeAnnotations import proteome_annotations
    from SchemaRefinery.SchemaAnnotation.AnnotationMerger import annotation_merger
    from SchemaRefinery.SchemaAnnotation.MatchSchemas import match_schemas
    from SchemaRefinery.utils.file_functions import check_and_make_directory


def main(args):

    check_and_make_directory(args.output_directory)

    # run annotation submodules
    if 'genbank' in args.annotation_options:
        genbank_file = genbank_annotations(args.input_files, args.schema_directory, args.output_directory, args.cpu_cores)

    if 'proteomes' in args.annotation_options:
        trembl, swiss = proteome_annotations(args.input_table, args.proteomes_directory, args.threads, args.retry, args.schema_directory, args.output_directory, args.cpu_cores)

    if 'matchSchemas' in args.annotation_options:
        matched_schemas = match_schemas(args.query_schema, args.subject_schema, args.output_directory, args.blast_score_ratio, args.cpu_cores)

    annotation_merger(args.uniprot_species, args.uniprot_genus, genbank_file, trembl, swiss, args.match_to_add, args.old_schema_columns, matched_schemas, args.output_directory)
