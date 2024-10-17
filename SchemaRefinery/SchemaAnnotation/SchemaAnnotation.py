import os
import pandas as pd
from argparse import Namespace
from functools import reduce
from typing import List, Optional, Tuple, Dict

try:
    from SchemaAnnotation import (proteome_fetcher as pf,
                                  proteome_splitter as ps,
                                  proteome_matcher as pm,
                                  genbank_annotations as ga,
                                  match_schemas as ms)
    from utils import (file_functions as ff)
except ModuleNotFoundError:
    from SchemaRefinery.SchemaAnnotation import (proteome_fetcher as pf,
                                                proteome_splitter as ps,
                                                proteome_matcher as pm,
                                                genbank_annotations as ga,
                                                match_schemas as ms)
    from SchemaRefinery.utils import (file_functions as ff)

def main(args: Namespace) -> None:
    """
    Main function to process proteome annotations based on provided arguments.

    Parameters
    ----------
    args : Namespace
        Parsed command-line arguments.

    Returns
    -------
    None
    """
    # Create the output directory if it doesn't exist
    ff.create_directory(args.output_directory)

    # Initialize a list to store result files
    results_files: List[str] = []

    # Add Chewie annotations to the results files if provided
    if args.chewie_annotations:
        results_files.extend(args.chewie_annotations)

    # Check if 'uniprot-proteomes' is in the annotation options
    if 'uniprot-proteomes' in args.annotation_options:
        uniprot_annotations_folder: str = os.path.join(args.output_directory, 'uniprot_annotations')
        # Fetch proteome data and store the directory path
        proteomes_directory: Optional[str] = pf.proteome_fetcher(args.proteome_table,
                                                                 uniprot_annotations_folder,
                                                                 args.threads,
                                                                 args.retry)

        if proteomes_directory is not None:
            # Split proteome records into TrEMBL and Swiss-Prot records
            split_data: Tuple[str, str, str, Dict[str, List[str]]] = ps.proteome_splitter(proteomes_directory,
                                                                    uniprot_annotations_folder)
            tr_file: str
            sp_file: str
            descriptions_file: str
            proteome_file_ids: Dict[str, List[str]]
            tr_file, sp_file, descriptions_file, proteome_file_ids = split_data

            # Align loci against proteome records
            annotations: List[str] = pm.proteome_matcher([tr_file, sp_file, descriptions_file],
                                                         proteome_file_ids,
                                                         args.schema_directory,
                                                         uniprot_annotations_folder,
                                                         args.cpu,
                                                         args.bsr,
                                                         args.translation_table,
                                                         args.clustering_sim,
                                                         args.clustering_cov,
                                                         args.size_ratio,
                                                         args.run_mode)
            results_files.extend(annotations)

    # Check if 'genbank' is in the annotation options
    if 'genbank' in args.annotation_options:
        genbank_annotation_folder: str = os.path.join(args.output_directory, 'genbank_annotations')
        # Process GenBank annotations
        genbank_file: str = ga.genbank_annotations(args.genbank_files,
                                                   args.schema_directory,
                                                   genbank_annotation_folder,
                                                   args.cpu,
                                                   args.bsr,
                                                   args.translation_table,
                                                   args.clustering_sim,
                                                   args.clustering_cov,
                                                   args.size_ratio,
                                                   args.run_mode)
        results_files.append(genbank_file)

    matched_schemas: Optional[str] = None
    # Check if 'match-schemas' is in the annotation options
    if 'match-schemas' in args.annotation_options:
        matched_schemas_folder: str = os.path.join(args.output_directory, 'matched_schemas')
        # Match schemas
        matched_schemas = ms.match_schemas(args.schema_directory,
                                           args.subject_schema,
                                           matched_schemas_folder,
                                           args.bsr,
                                           args.cpu)
        results_files.append(matched_schemas)

    # Merge all results into a single file
    dfs: List[pd.DataFrame] = []
    for file in results_files:
        current_df: pd.DataFrame = pd.read_csv(file, delimiter='\t', dtype=str)
        if 'Locus' not in current_df.columns:
            current_df = current_df.rename({'Locus_ID': 'Locus'}, axis=1)
        dfs.append(current_df)

    if args.subject_annotations and matched_schemas:
        # Read TSV with subject schema annotations
        if args.subject_columns:
            columns: List[str] = ["Locus" if col == "Locus_ID" else col for col in args.subject_columns]
            if "Locus" not in columns:
                columns = ['Locus'] + columns

            match_add: pd.DataFrame = pd.read_csv(args.subject_annotations, delimiter='\t',
                                                  usecols=columns, dtype=str)
        else:
            match_add: pd.DataFrame = pd.read_csv(args.subject_annotations, delimiter='\t', dtype=str)

        # Merge columns so that both table to add and reference have locus_ID
        merged_match: pd.DataFrame = pd.merge(match_add, dfs[-1], on='Locus',
                                              how='left').fillna('')
        merged_match = merged_match[merged_match.columns.tolist()[1:]]

        dfs[-1] = merged_match

    # Merge all dataframes based on locus identifier
    merged_table: pd.DataFrame = reduce(lambda a, b: pd.merge(a, b, on=['Locus'],
                                                              how='left'), dfs).fillna('')

    # Save the merged table to a TSV file
    merged_table.to_csv(os.path.join(args.output_directory, 'merged_file.tsv'),
                        sep='\t', index=False)
