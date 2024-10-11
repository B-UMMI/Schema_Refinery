#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Purpose
-------

This module enables the whole process of creating the
Genbank, TrEMBL and Swiss-Prot annotations and merging them
alongside the UniProtFinder and/or match_schemas annotations into
a single TSV file.

Code documentation
------------------
"""


import os
from functools import reduce

import pandas as pd

try:
    from SchemaAnnotation import (proteome_fetcher as pf,
                                  proteome_splitter as ps,
                                  proteome_matcher as pm,
                                  genbank_annotations as ga,
                                  match_schemas as ms)
    from utils import file_functions as ff
except ModuleNotFoundError:
    from SchemaRefinery.SchemaAnnotation import (proteome_fetcher as pf,
                                                proteome_splitter as ps,
                                                proteome_matcher as pm,
                                                genbank_annotations as ga,
                                                match_schemas as ms)
    from SchemaRefinery.utils import file_functions as ff


def main(args):

    ff.create_directory(args.output_directory)

    results_files = []

    if args.chewie_annotations:
        results_files.extend(args.chewie_annotations)

    if 'uniprot-proteomes' in args.annotation_options:
        proteomes_directory = pf.proteome_fetcher(args.proteome_table,
                                                  args.output_directory,
                                                  args.threads,
                                                  args.retry)

        if proteomes_directory is not None:
            # Split proteome records into two files, one with TrEMBL records
            # and another with Swiss-Prot records
            split_data = ps.proteome_splitter(proteomes_directory,
                                              args.output_directory)
            tr_file, sp_file, descriptions_file = split_data

            # Align loci against proteome records
            annotations = pm.proteome_matcher(args.schema_directory,
                                              [tr_file, sp_file, descriptions_file],
                                              args.output_directory,
                                              args.cpu,
                                              args.bsr)
            results_files.extend(annotations)
    if 'genbank' in args.annotation_options:
        genbank_file = ga.genbank_annotations(args.genbank_files,
                                              args.schema_directory,
                                              args.output_directory,
                                              args.cpu,
                                              args.bsr,
                                              args.translation_table,
                                              args.clustering_sim,
                                              args.clustering_cov,
                                              args.size_ratio,
                                              args.run_mode)
        results_files.append(genbank_file)
    matched_schemas = None
    if 'match-schemas' in args.annotation_options:
        matched_schemas = ms.match_schemas(args.schema_directory,
                                           args.subject_schema,
                                           args.output_directory,
                                           args.bsr,
                                           args.cpu)
        results_files.append(matched_schemas)

    # Merge all results into a single file
    dfs = []
    for file in results_files:
        current_df = pd.read_csv(file, delimiter='\t', dtype=str)
        if 'Locus' not in current_df.columns:
            current_df = current_df.rename({'Locus_ID': 'Locus'}, axis=1)
        dfs.append(current_df)

    if args.subject_annotations and matched_schemas:
        # Read TSV with subject schema annotations
        if args.subject_columns:
            columns = ["Locus" if col == "Locus_ID" else col for col in args.subject_columns]
            if "Locus" not in columns:
                columns = ['Locus'] + columns

            match_add = pd.read_csv(args.subject_annotations, delimiter='\t',
                                    usecols=columns, dtype=str)
        else:
            match_add = pd.read_csv(args.subject_annotations, delimiter='\t', dtype=str)

        # Merge columns so that both table to add and reference have locus_ID
        merged_match = pd.merge(match_add, dfs[-1], on='Locus',
                                how='left').fillna('')
        merged_match = merged_match[merged_match.columns.tolist()[1:]]

        dfs[-1] = merged_match

    # Merge all dataframes based on locus identifier
    merged_table = reduce(lambda a, b: pd.merge(a, b, on=['Locus'],
                                                how='left'), dfs).fillna('')

    merged_table.to_csv(os.path.join(args.output_directory, 'merged_file.tsv'),
                        sep='\t', index=False)
