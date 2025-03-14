import pandas as pd
from typing import Dict, Tuple, List

try:
    from utils import (pandas_functions as upf,
					   print_functions as prf)
except ModuleNotFoundError:
    from SchemaRefinery.utils import (pandas_functions as upf,
									  print_functions as prf)

def consolidate_annotations(consolidate_annotation_files: List[str],
                            cleanup: bool,
                            output_file: str) -> str:
    """
    Consolidates annotations from different annotation files into a single file.
    
    Parameters
    ----------
    consolidate_annotation_files: List[str]:
        List of files of annotations to consolidate.
    cleanup: bool
        If the final file will or not have duplicates. Advised for the use of match schemas annotations.
    output_file: str
        Path to the directory where output files will be saved.

    Returns
    -------
    str
        Path to the annotations file.
    """

    
    first_df = pd.read_csv(consolidate_annotation_files[0], delimiter='\t', dtype=str, index_col=False)
    second_df = pd.read_csv(consolidate_annotation_files[1], delimiter='\t', dtype=str, index_col=False)
    
    first_filtered = first_df[(first_df.iloc[:, 0] != 'Not Matched') & (first_df.iloc[:, 1] != 'Not Matched')]
    second_filtered = second_df[(second_df.iloc[:, 0] != 'Not Matched') & (second_df.iloc[:, 1] != 'Not Matched')]

    # f0=s0
    if first_filtered.iloc[:, 0].equals(second_filtered.iloc[:, 0]):
        upf.merge_files_by_column_values(consolidate_annotation_files[0],
                                            consolidate_annotation_files[1],
                                            0,
                                            0,
                                            output_file)

    # f0=s1
    elif first_filtered.iloc[:, 0].equals(second_filtered.iloc[:, 1]):
        upf.merge_files_by_column_values(consolidate_annotation_files[0],
                                            consolidate_annotation_files[1],
                                            0,
                                            1,
                                            output_file)

    # f1=s0
    elif first_filtered.iloc[:, 1].equals(second_filtered.iloc[:, 0]):
        upf.merge_files_by_column_values(consolidate_annotation_files[0],
                                            consolidate_annotation_files[1],
                                            1,
                                            0,
                                            output_file)

    # f1=s1
    elif first_filtered.iloc[:, 1].equals(second_filtered.iloc[:, 1]):
        upf.merge_files_by_column_values(consolidate_annotation_files[0],
                                            consolidate_annotation_files[1],
                                            1,
                                            1,
                                            output_file)



    ## depois com o ficheiro 2 começa o loop atá ao fim i < len(list) e i > 1
    ## loop com output_file como input e list[i] como segundo input
    for i in range(2, len(consolidate_annotation_files)):
        old_df = pd.read_csv(output_file, delimiter='\t', dtype=str, index_col=False)
        new_df = pd.read_csv(consolidate_annotation_files[i], delimiter='\t', dtype=str, index_col=False)
        
        old_filtered = old_df[(old_df.iloc[:, 0] != 'Not Matched') & (old_df.iloc[:, 1] != 'Not Matched')]
        new_filtered = new_df[(new_df.iloc[:, 0] != 'Not Matched') & (new_df.iloc[:, 1] != 'Not Matched')]

        # f0=s0
        if old_filtered.iloc[:, 0].equals(new_filtered.iloc[:, 0]):
            upf.merge_files_by_column_values(output_file,
                                                consolidate_annotation_files[i],
                                                0,
                                                0,
                                                output_file)

        # f0=s1
        elif old_filtered.iloc[:, 0].equals(new_filtered.iloc[:, 1]):
            upf.merge_files_by_column_values(output_file,
                                                consolidate_annotation_files[i],
                                                0,
                                                1,
                                                output_file)

        # f1=s0
        elif old_filtered.iloc[:, 1].equals(new_filtered.iloc[:, 0]):
            upf.merge_files_by_column_values(output_file,
                                                consolidate_annotation_files[i],
                                                1,
                                                0,
                                                output_file)

        # f1=s1
        elif old_filtered.iloc[:, 1].equals(new_filtered.iloc[:, 1]):
            upf.merge_files_by_column_values(output_file,
                                                consolidate_annotation_files[i],
                                                1,
                                                1,
                                                output_file)

    # REVER
    # If cleanup TRUE dedeup the final file
    # Right now only works for ms
    if cleanup:
        ann_df = pd.read_csv(output_file, delimiter='\t', dtype=str, index_col=False)
        ann_df.sort_values(by=ann_df.columns[0], ascending=False, inplace=True)
        ann_df.drop_duplicates(subset=ann_df.columns[0], inplace=True)
        ann_df.to_csv(output_file, sep="\t", index=False)

    return output_file