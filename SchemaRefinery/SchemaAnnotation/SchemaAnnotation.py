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

    # Initialize a list to store result files
    results_files: List[str] = []

    # Check if 'uniprot-proteomes' is in the annotation options
    if 'uniprot-proteomes' in args.annotation_options:
        prf.print_message('Running Annotation with proteomes.', 'info')
        uniprot_annotations_folder: str = os.path.join(args.output_directory, 'uniprot_annotations')
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
        matched_annotations = os.path.join(args.output_directory, "matched_annotations.tsv")

        matched_df = pd.read_csv(args.matched_schemas, delimiter='\t', dtype=str, index_col=False)
        annotations_df = pd.read_csv(args.subject_annotations, delimiter='\t', dtype=str, index_col=False)

        matched_0_filtered = matched_df[matched_df.iloc[:, 0] != 'Not matched'].sort_values(by=matched_df.columns[0]).drop_duplicates(subset=['Query']).reset_index(drop=True)
        matched_1_filtered = matched_df[matched_df.iloc[:, 1] != 'Not matched'].sort_values(by=matched_df.columns[1]).drop_duplicates(subset=['Subject']).reset_index(drop=True)
        annotations_sorted = annotations_df.sort_values(by=annotations_df.columns[0]).reset_index(drop=True)

        prf.print_message(f'{matched_0_filtered.iloc[:, 0]}', 'info')
        prf.print_message(f'{annotations_sorted.iloc[:, 0]}', 'info')

        if matched_0_filtered.iloc[:, 0].equals(annotations_sorted.iloc[:, 0]):
            prf.print_message("Annotating from the Query", "info")
            upf.merge_files_by_column_values(args.matched_schemas,
                                            args.subject_annotations,
                                            0,
                                            0,
                                            matched_annotations)
        if matched_1_filtered.iloc[:, 1].equals(annotations_df.iloc[:, 0]):
            prf.print_message('Annotating from the Subject', 'info')                                    
            upf.merge_files_by_column_values(args.matched_schemas,
                                            args.subject_annotations,
                                            1,
                                            0,
                                            matched_annotations)
        else:
            prf.print_message('No matches found in columns', 'info')
            
        results_files.append(matched_annotations)
        prf.print_message('Matching successfully completed.', 'info')

    if 'consolidate' in args.annotation_options:
        prf.print_message("Consolidating annoations...", "info")
        consolidated_annotations = os.path.join(args.output_directory, "consolidated_annotations.tsv")
        consolidated_annotations_final: str = cs.consolidate_annotations(args.consolidate_annotations,
                                                                            args.consolidate_cleanup,
                                                                            consolidated_annotations)

        results_files.append(consolidated_annotations_final)
        prf.print_message('Annotation consolidation successfully completed.', 'info')

    # Add Chewie annotations to the results files if provided
    if args.chewie_annotations:
        results_files.extend(args.chewie_annotations)

    merged_file_path = os.path.join(args.output_directory, 'annotations_summary.tsv')

    # If only one result file is present, copy it to the output directory
    if len(results_files) == 1:
        shutil.copy(results_files[0], merged_file_path)
    else:
        # Merge all results into a single file
        upf.merge_files_into_same_file_by_key(results_files, 'Locus', merged_file_path)

    if args.best_annotations_bsr:
        priority_dict = {}
        if 'genbank' in args.annotation_options:
            priority_dict.update({
                'Genbank_BSR' : ['Locus', 'Genbank_ID', 'Genbank_gene_name', 'Genbank_product', 'Genbank_BSR']
            })
        if 'uniprot-proteomes' in args.annotation_options:
            priority_dict.update({
                'Proteome_BSR' : ['Locus', 'Proteome_ID', 'Proteome_product', 'Proteome_gene_name', 'Proteome_BSR']
            })
        # Process the merged file based on the priority dictionary
        output_file = os.path.join(args.output_directory, 'best_annotations_user_input.tsv')
        # Define the columns to include in the output file
        output_columns = ['Locus', 'Protein_ID', 'Protein_gene_name', 'Protein_product', 'Protein_BSR', 'Source']
        upf.process_tsv_with_priority(merged_file_path, priority_dict, output_file, args.best_annotations_bsr, output_columns)
    # Clean up temporary files
    if not args.no_cleanup:
        if 'match-schemas' not in args.annotation_options:
            prf.print_message("Cleaning up temporary files...", 'info')
        # Remove temporary files
        ff.cleanup(args.output_directory, [merged_file_path, output_file, logf.get_log_file_path(gb.LOGGER)])