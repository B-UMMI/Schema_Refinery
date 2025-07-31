import os
import shutil
import pandas as pd
from argparse import Namespace
from typing import List, Optional, Tuple, Dict


try:
	from SchemaAnnotation import (proteome_fetcher as pf,
								  proteome_splitter as ps,
								  proteome_matcher as pm,
								  genbank_annotations as ga,
                                  consolidate as cs)
	from utils import (file_functions as ff,
					   pandas_functions as upf,
					   print_functions as prf,
					   logger_functions as logf,
					   globals as gb)
except ModuleNotFoundError:
	from SchemaRefinery.SchemaAnnotation import (proteome_fetcher as pf,
												proteome_splitter as ps,
												proteome_matcher as pm,
												genbank_annotations as ga,
                                                consolidate as cs)
	from SchemaRefinery.utils import (file_functions as ff,
									  pandas_functions as upf,
									  print_functions as prf,
									  logger_functions as logf,
									  globals as gb)


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
    output_d= os.path.abspath(args.output_directory)

    # Initialize a list to store result files
    results_files: List[str] = []

    # Check if 'uniprot-proteomes' is in the annotation options
    if 'uniprot-proteomes' in args.annotation_options:
        prf.print_message('Running Annotation with proteomes.', 'info')
        uniprot_annotations_folder: str = os.path.join(output_d, 'uniprot_annotations')
        merged_file_path = os.path.join(output_d, 'uniprot_annotations.tsv')
        # Fetch proteome data and store the directory path
        proteomes_directory: Optional[str] = pf.proteome_fetcher(args.proteome_table,
                                                                 uniprot_annotations_folder,
                                                                 args.threads,
                                                                 args.retry)

        if proteomes_directory is not None:
            # Split proteome records into TrEMBL and Swiss-Prot records
            prf.print_message('Spliting the annotations into Swiss and TrEMBL.', 'info')
            split_data: Tuple[str, str, str, Dict[str, List[str]]] = ps.proteome_splitter(proteomes_directory,
                                                                    uniprot_annotations_folder)
            tr_file: str
            sp_file: str
            descriptions_file: str
            proteome_file_ids: Dict[str, List[str]]
            tr_file, sp_file, descriptions_file, proteome_file_ids = split_data

            # Align loci against proteome records
            prf.print_message('Matching the annotations.', 'info')
            annotations: List[str] = pm.proteome_matcher([tr_file, sp_file, descriptions_file],
                                                         proteome_file_ids,
                                                         args.schema_directory,
                                                         uniprot_annotations_folder,
                                                         args.cpu,
                                                         args.bsr,
                                                         args.translation_table,
                                                         args.run_mode,
                                                         args.proteome_ids_to_add)
            results_files.extend(annotations)
            prf.print_message('Matching successfully completed.', 'info')
    

    # Check if 'genbank' is in the annotation options
    if 'genbank' in args.annotation_options:
        prf.print_message('Running Annotation with GenBank.', 'info')
        genbank_annotation_folder: str = os.path.join(output_d, 'genbank_annotations')
        merged_file_path = os.path.join(output_d, 'genbank_annotations.tsv')
        # Process GenBank annotations
        genbank_file: str = ga.genbank_annotations(args.genbank_files,
                                                   args.schema_directory,
                                                   genbank_annotation_folder,
                                                   args.cpu,
                                                   args.bsr,
                                                   args.translation_table,
                                                   args.run_mode,
                                                   args.extra_genbank_table_columns,
                                                   args.genbank_ids_to_add)
        results_files.append(genbank_file)
        prf.print_message('Matching successfully completed.', 'info')

    matched_schemas: Optional[str] = None
    # Check if 'match-schemas' is in the annotation options
    if 'match-schemas' in args.annotation_options:
        # Merge matched loci with their annotation
        prf.print_message('Running Annotation with MatchSchemas.', 'info')
        merged_file_path = os.path.join(output_d, "matched_annotations.tsv")
        matched_annotations = None

        # Create df from the tsv files given and compare the sorted and filtered columns
        matched_df = pd.read_csv(args.matched_schemas, delimiter='\t', dtype=str, index_col=False)
        annotations_df = pd.read_csv(args.match_annotations, delimiter='\t', dtype=str, index_col=False)

        matched_0_filtered = matched_df[matched_df.iloc[:, 0] != 'Not matched'].sort_values(by=matched_df.columns[0]).drop_duplicates(subset=['Query']).reset_index(drop=True)
        matched_1_filtered = matched_df[matched_df.iloc[:, 1] != 'Not matched'].sort_values(by=matched_df.columns[1]).drop_duplicates(subset=['Subject']).reset_index(drop=True)
        annotations_sorted = annotations_df.sort_values(by=annotations_df.columns[0]).reset_index(drop=True)

        matches = {
            'f0s0': matched_0_filtered.iloc[:, 0].isin(annotations_sorted.iloc[:, 0]).sum(),
            'f1s0': matched_1_filtered.iloc[:, 1].isin(annotations_sorted.iloc[:, 0]).sum(),
            }

        best_match = max(matches, key=matches.get)

        # Depending on which columns are a match run different versions of the merging
        if best_match == 'f0s0':
            prf.print_message("Annotating from the Query", "info")
            mismatched_f0s0 = matched_0_filtered.iloc[:, 0][~matched_0_filtered.iloc[:, 0].isin(annotations_sorted.iloc[:, 0])]
            prf.print_message("Mismatched rows in Query compared and annotations:")
            prf.print_message(f"From Query: {mismatched_f0s0}", 'info')
            mismatched_s0f0 = annotations_sorted.iloc[:, 0][~annotations_sorted.iloc[:, 0].isin(matched_0_filtered.iloc[:, 0])]
            prf.print_message(f"From Annotation: {mismatched_s0f0}", 'info')

            matched_annotations: str = upf.merge_files_by_column_values(args.matched_schemas,
                                            args.match_annotations,
                                            0,
                                            0,
                                            merged_file_path)
        elif best_match == 'f1s0':
            prf.print_message('Annotating from the Subject', 'info') 
            mismatched_f1s0 = matched_1_filtered.iloc[:, 1][~matched_1_filtered.iloc[:, 1].isin(annotations_sorted.iloc[:, 0])]
            prf.print_message("Mismatched rows in Subject compared and annotations:")
            prf.print_message(f"From Subject: {mismatched_f1s0}", 'info')
            mismatched_s0f1 = annotations_sorted.iloc[:, 0][~annotations_sorted.iloc[:, 0].isin(matched_1_filtered.iloc[:, 1])]
            prf.print_message(f"From Annotation: {mismatched_s0f1}", 'info')

            matched_annotations: str = upf.merge_files_by_column_values(args.matched_schemas,
                                            args.match_annotations,
                                            1,
                                            0,
                                            merged_file_path)
        else:
            prf.print_message('No matches found in columns', 'info')
            
        results_files.append(matched_annotations)
        prf.print_message('Matching successfully completed.', 'info')

    # Check if 'consolidate' is in the annotation options
    if 'consolidate' in args.annotation_options:
        prf.print_message("Consolidating annotations...", "info")
        merged_file_path = os.path.join(output_d, "consolidated_annotations.tsv")
        consolidated_annotations = None
        consolidated_annotations: str = cs.consolidate_annotations(args.consolidate_annotations,
                                    args.consolidate_cleanup,
                                    merged_file_path)

        results_files.append(consolidated_annotations)
        prf.print_message('Annotation consolidation successfully completed.', 'info')
        prf.print_message('')

    # Add Chewie annotations to the results files if provided
    if args.chewie_annotations:
        results_files.append(args.chewie_annotations)
        upf.merge_files_into_same_file_by_key(results_files, 'Locus', merged_file_path)

    # If only one result file is present, copy it to the output directory
    if len(results_files) == 1:
        if 'match-schemas' not in args.annotation_options and 'consolidate' not in args.annotation_options:
            shutil.copy(results_files[0], merged_file_path)
    else:
        # Merge all results into a single file
        upf.merge_files_into_same_file_by_key(results_files, 'Locus', merged_file_path)

    # Print final statistics
    annotations_count = 0
    hypoteticals = 0
    with open(merged_file_path, "r") as f:
        for line in f:
            parts = line.strip().split("\t")
            # Only count rows that are not empty
            # Match schemas mode always has 2 extra written rows even if there are no annotations
            if len(parts) > 3:
                # Only count if there no missing information
                if any(field != "NA" for field in parts[1:]):
                    annotations_count += 1
            # Count hypothetical proteins
            if any("hypothetical protein" in field.lower() for field in parts):
                hypoteticals += 1

    prf.print_message('')
    prf.print_message(f'A total of {annotations_count-1} loci were annotated.', 'info')
    prf.print_message(f'From these {hypoteticals} loci where annotated as "hypothetical proteins".', 'info')
    prf.print_message('')

    if not args.no_cleanup:
        if 'match-schemas' not in args.annotation_options:
            prf.print_message("Cleaning up temporary files...", 'info')
        ff.cleanup(output_d, [merged_file_path, logf.get_log_file_path(gb.LOGGER)])