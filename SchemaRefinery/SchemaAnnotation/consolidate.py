import pandas as pd
from typing import Dict, Tuple, List

try:
    from utils import (pandas_functions as upf,
					   print_functions as prf)
except ModuleNotFoundError:
    from SchemaRefinery.utils import (pandas_functions as upf,
									  print_functions as prf)

def consolidate_annotations(uniprot_annotations_file: str,
                            genbank_annotations_file: str,
                            ms_annotations_file: str,
                            cleanup: bool,
                            output_file: str) -> str:
    """
    Consolidates annotations from different annotation files into a single file.
    
    Parameters
    ----------
    uniprot_annotations_file: str
        Path to the uniprot annotations file.
    genbank_annotations_file: str
        Path to the genbank annotations file.
    ms_annotations_file: str
        Path to the match schemas annotations file.
    cleanup: bool
        If the final file will or not have duplicates. Advised for the use of match schemas annotations.
    output_file: str
        Path to the directory where output files will be saved.

    Returns
    -------
    str
        Path to the annotations file.
    """

    # Based on which files are not none?? --> Yup
    # Merge uniprot and genebank
    if uniprot_annotations_file and genbank_annotations_file:
        if ms_annotations_file is None:
            ## use merge file based on column
            upf.merge_files_by_column_values(uniprot_annotations_file, 
                                                genbank_annotations_file, 
                                                0, 
                                                0, 
                                                output_file)
        continue
    
    # Merge genbank and ms  
    if ms_annotations_file and genbank_annotations_file:
        if uniprot_annotations_file is None:
            ## compare locus column with subject and query (no Not Matched)
            ms_df = pd.read_csv(ms_annotations_file, delimiter='\t', dtype=str, index_col=False)
	        genbank_df = pd.read_csv(genbank_annotations_file, delimiter='\t', dtype=str, index_col=False)
            ms_query_filtered = ms_df[ms_df['Query'] != 'Not Matched']
            ms_subject_filtered = ms_df[ms_df['Subject'] != 'Not Matched']
            if ms_query_filtered['Query'].equals(genbank_df['Locus']):
                upf.merge_files_by_column_values(ms_annotations_file,
                                                    genbank_annotations_file,
                                                    0,
                                                    0,
                                                    output_file)
            if ms_subject_filtered['Subject'].equals(genbank_df['Locus']):
                upf.merge_files_by_column_values(ms_annotations_file,
                                                    genbank_annotations_file,
                                                    0,
                                                    0,
                                                    output_file)
        continue
    
    # Merge uniprot and ms
    if ms_annotations_file and uniprot_annotations_file:
        if genbank_annotations_file is None:
            ## compare locus column with subject and query (no Not Matched)
            ms_df = pd.read_csv(ms_annotations_file, delimiter='\t', dtype=str, index_col=False)
	        uniprot_df = pd.read_csv(uniprot_annotations_file, delimiter='\t', dtype=str, index_col=False)
            ms_query_filtered = ms_df[ms_df['Query'] != 'Not Matched']
            ms_subject_filtered = ms_df[ms_df['Subject'] != 'Not Matched']
            if ms_query_filtered['Query'].equals(uniprot_df['Locus']):
                upf.merge_files_by_column_values(ms_annotations_file,
                                                    uniprot_annotations_file,
                                                    0,
                                                    0,
                                                    output_file)
            if ms_subject_filtered['Subject'].equals(uniprot_df['Locus']):
                upf.merge_files_by_column_values(ms_annotations_file,
                                                    uniprot_annotations_file,
                                                    1,
                                                    0,
                                                    output_file)
        continue

    # Merge uniprot and genbank and then MS
    if ms_annotations_file and uniprot_annotations_file and genbank_annotations_file:
        upf.merge_files_by_column_values(uniprot_annotations_file, 
                                                genbank_annotations_file, 
                                                0, 
                                                0, 
                                                unigen_annotations_file)
        ## compare locus column with subject and query (no Not Matched)
        ms_df = pd.read_csv(ms_annotations_file, delimiter='\t', dtype=str, index_col=False)
        unigen_df = pd.read_csv(unigen_annotations_file, delimiter='\t', dtype=str, index_col=False)
        ms_query_filtered = ms_df[ms_df['Query'] != 'Not Matched']
        ms_subject_filtered = ms_df[ms_df['Subject'] != 'Not Matched']
        if ms_query_filtered['Query'].equals(unigen_df['Locus']):
            upf.merge_files_by_column_values(ms_annotations_file,
                                                unigen_annotations_file,
                                                0,
                                                0,
                                                consolidated_annotation_file)
        if ms_subject_filtered['Subject'].equals(unigen_df['Locus']):
            upf.merge_files_by_column_values(ms_annotations_file,
                                                unigen_annotations_file,
                                                1,
                                                0,
                                                output_file)
    continue


    # If cleanup TRUE dedeup the final file
    # Right now only works for ms
    if cleanup:
        ann_df = pd.read_csv(output_file, delimiter='\t', dtype=str, index_col=False)
        ann_df.sort_values(by='Query', ascending=False)
        ann_df.drop_duplicates(subset=['Query'], inplace=TRUE)
        ann_df.to_csv(output_file, sep="\t")

    return output_file