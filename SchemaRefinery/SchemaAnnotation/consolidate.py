import os
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

    # Create dataframes out of the first 2 files in the list
    prf.print_message('Format tsv files for merging...', 'info')

    with open(output_file, 'w') as f:
        f.write("")

    # Check if all files exist
    for i in range(len(consolidate_annotation_files)):
        if os.path.getsize(consolidate_annotation_files[i]) == 0:
            prf.print_message(f'The file number {i+1} is empty.', 'error')
            sys.exit()
        else:
            prf.print_message(f'File {i+1}: {consolidate_annotation_files[i]}', 'info')
    
    prf.print_message('')
    
    first_df = pd.read_csv(consolidate_annotation_files[0], delimiter='\t', dtype=str, index_col=False)
    second_df = pd.read_csv(consolidate_annotation_files[1], delimiter='\t', dtype=str, index_col=False)
    
    # Filter and sort the first and second column of both files
    first_0_filtered = first_df[first_df.iloc[:, 0] != 'Not matched'].sort_values(by=first_df.columns[0]).drop_duplicates(subset=first_df.columns[0]).reset_index(drop=True)
    second_0_filtered = second_df[(second_df.iloc[:, 0] != 'Not matched')].sort_values(by=second_df.columns[0]).drop_duplicates(subset=second_df.columns[0]).dropna().reset_index(drop=True)
    first_1_filtered = first_df[(first_df.iloc[:, 1] != 'Not matched')].sort_values(by=first_df.columns[1]).drop_duplicates(subset=first_df.columns[1]).reset_index(drop=True)
    second_1_filtered = second_df[(second_df.iloc[:, 1] != 'Not matched')].sort_values(by=second_df.columns[1]).drop_duplicates(subset=second_df.columns[1]).dropna().reset_index(drop=True)

    matches = {
    'f0s0': first_0_filtered.iloc[:, 0].isin(second_0_filtered.iloc[:, 0]).sum(),
    'f0s1': first_0_filtered.iloc[:, 0].isin(second_1_filtered.iloc[:, 1]).sum(),
    'f1s0': first_1_filtered.iloc[:, 1].isin(second_0_filtered.iloc[:, 0]).sum(),
    'f1s1': first_1_filtered.iloc[:, 1].isin(second_1_filtered.iloc[:, 1]).sum()
    }

    best_match = max(matches, key=matches.get)

    if matches[best_match] == 0:
        prf.print_message('No matches were found. The final file will have no annotations.', 'warning')

    prf.print_message('Merging first 2 files...', 'info')
    # Compare the first columns of both files
    if best_match == 'f0s0':
        prf.print_message('First columns are a match', 'info')
        mismatched_f0s0 = first_0_filtered.iloc[:, 0][~first_0_filtered.iloc[:, 0].isin(second_0_filtered.iloc[:, 0])]
        prf.print_message(f"Ids from the First file that didn't match: {mismatched_f0s0}", 'info')
        mismatched_s0f0 = second_0_filtered.iloc[:, 0][~second_0_filtered.iloc[:, 0].isin(first_0_filtered.iloc[:, 0])]
        prf.print_message(f"Ids from the second file that didn't match: {mismatched_s0f0}", 'info')
        upf.merge_files_by_column_values_df(first_df,
                                            second_0_filtered,
                                            0,
                                            0,
                                            output_file,
                                            '_file_1',
                                            '_file_2')
    # Compare the first column of the first file and the second one of the second file
    elif best_match == 'f0s1':
        prf.print_message('First and second columns are a match', 'info')
        mismatched_f0s1 = first_0_filtered.iloc[:, 0][~first_0_filtered.iloc[:, 0].isin(second_1_filtered.iloc[:, 1])]
        prf.print_message(f"Ids from the First file that didn't match: {mismatched_f0s1}", 'info')
        mismatched_s1f0 = second_1_filtered.iloc[:, 1][~second_1_filtered.iloc[:, 1].isin(first_0_filtered.iloc[:, 0])]
        prf.print_message(f"Ids from the Second file that didn't match: {mismatched_s1f0}", 'info')
        upf.merge_files_by_column_values_df(second_df,
                                            first_0_filtered,
                                            1,
                                            0,
                                            output_file,
                                            '_file_2',
                                            '_file_1')
        # Make the first column the one that is common among all 
        out_df = pd.read_csv(output_file, delimiter='\t', dtype=str, index_col=False)
        cols = list(out_df.columns)
        out_df = out_df[[cols[1], cols[0]] + cols[2:]]
        out_df.to_csv(output_file, sep='\t', index=False)
    # Compare the second column of the first file and the first one of the second file
    elif best_match == 'f1s0':
        prf.print_message('Second and first columns are a match', 'info')
        mismatched_f1s0 = first_1_filtered.iloc[:, 1][~first_1_filtered.iloc[:, 1].isin(second_0_filtered.iloc[:, 0])]
        prf.print_message(f"Ids from the First file that didn't match: {mismatched_f1s0}", 'info')
        mismatched_s0f0 = second_0_filtered.iloc[:, 0][~second_0_filtered.iloc[:, 0].isin(first_0_filtered.iloc[:, 0])]
        prf.print_message(f"Ids from the Second file that didn't match: {mismatched_s0f0}", 'info')
        upf.merge_files_by_column_values_df(first_df,
                                            second_0_filtered,
                                            1,
                                            0,
                                            output_file,
                                            '_file_1',
                                            '_file_2')
        # Make the first column the one that is common among all 
        out_df = pd.read_csv(output_file, delimiter='\t', dtype=str, index_col=False)
        cols = list(out_df.columns)
        out_df = out_df[[cols[1], cols[0]] + cols[2:]]
        out_df.to_csv(output_file, sep='\t', index=False)
    # Compare the second columns of both files
    elif best_match == 'f1s1':
        prf.print_message('Second columns are a match', 'info')
        mismatched_f1s1 = first_1_filtered.iloc[:, 1][~first_1_filtered.iloc[:, 1].isin(second_1_filtered.iloc[:, 1])]
        prf.print_message(f"Ids from the First file that didn't match: {mismatched_f1s1}", 'info')
        mismatched_s1f1 = second_1_filtered.iloc[:, 1][~second_1_filtered.iloc[:, 1].isin(first_1_filtered.iloc[:, 1])]
        prf.print_message(f"Ids from the Second file that didn't match: {mismatched_s1f1}", 'info')
        upf.merge_files_by_column_values_df(first_df,
                                            second_1_filtered,
                                            1,
                                            1,
                                            output_file,
                                            '_file_1',
                                            '_file_2')
        # Make the first column the one that is common among all 
        out_df = pd.read_csv(output_file, delimiter='\t', dtype=str, index_col=False)
        cols = list(out_df.columns)
        out_df = out_df[[cols[1], cols[0]] + cols[2:]]
        out_df.to_csv(output_file, sep='\t', index=False)
    else:
        prf.print_message('No columns are a match', 'info')

    # Merge the remaining files
    if len(consolidate_annotation_files) > 2:
        prf.print_message('Merging the remaining files...', 'info')
        for i in range(2, len(consolidate_annotation_files)):
            prf.print_message(f"Merging file: {i+1}/{len(consolidate_annotation_files)}", "info", end='\r', flush=True)
            old_df = pd.read_csv(output_file, delimiter='\t', dtype=str, index_col=False)
            new_df = pd.read_csv(consolidate_annotation_files[i], delimiter='\t', dtype=str, index_col=False)
            
            # Filter and sort the first and second column of both files
            old_0_filtered = old_df[(old_df.iloc[:, 0] != 'Not matched')].sort_values(by=old_df.columns[0]).drop_duplicates(subset=old_df.columns[0]).reset_index(drop=True)
            new_0_filtered = new_df[(new_df.iloc[:, 0] != 'Not matched')].sort_values(by=new_df.columns[0]).drop_duplicates(subset=new_df.columns[0]).reset_index(drop=True)
            old_1_filtered = old_df[(old_df.iloc[:, 1] != 'Not matched')].sort_values(by=old_df.columns[1]).drop_duplicates(subset=old_df.columns[1]).reset_index(drop=True)
            new_1_filtered = new_df[(new_df.iloc[:, 1] != 'Not matched')].sort_values(by=new_df.columns[1]).drop_duplicates(subset=new_df.columns[1]).reset_index(drop=True)

            matches = {
            'f0s0': old_0_filtered.iloc[:, 0].isin(new_0_filtered.iloc[:, 0]).sum(),
            'f0s1': old_0_filtered.iloc[:, 0].isin(new_1_filtered.iloc[:, 1]).sum(),
            'f1s0': old_1_filtered.iloc[:, 1].isin(new_0_filtered.iloc[:, 0]).sum(),
            'f1s1': old_1_filtered.iloc[:, 1].isin(new_1_filtered.iloc[:, 1]).sum()
            }

            best_match = max(matches, key=matches.get)


            # Compare the first columns of both files
            if best_match == 'f0s0':
                prf.print_message('First columns are a match', 'info')
                mismatched_f0s0 = old_0_filtered.iloc[:, 0][~old_0_filtered.iloc[:, 0].isin(new_0_filtered.iloc[:, 0])]
                prf.print_message(f"Ids from the consolidated file that didn't match: {mismatched_f0s0}", 'info')
                mismatched_s0f0 = new_0_filtered.iloc[:, 0][~new_0_filtered.iloc[:, 0].isin(old_0_filtered.iloc[:, 0])]
                prf.print_message(f"Ids from the {i+1} file that didn't match: {mismatched_s0f0}", 'info')
                upf.merge_files_by_column_values_df(old_df,
                                                    new_0_filtered,
                                                    0,
                                                    0,
                                                    output_file,
                                                    None,
                                                    f'_file_{i+1}')
            # Compare the first column of the first file and the second one of the second file
            elif best_match == 'f0s1':
                prf.print_message('First and second columns are a match', 'info')
                mismatched_f0s1 = old_0_filtered.iloc[:, 0][~old_0_filtered.iloc[:, 0].isin(new_1_filtered.iloc[:, 1])]
                prf.print_message(f"Ids from the consolidated file that didn't match: {mismatched_f0s1}", 'info')
                mismatched_s1f0 = new_1_filtered.iloc[:, 1][~new_1_filtered.iloc[:, 1].isin(old_0_filtered.iloc[:, 0])]
                prf.print_message(f"Ids from the {i+1} file that didn't match: {mismatched_s1f0}", 'info')
                upf.merge_files_by_column_values_df(new_df,
                                                    old_0_filtered,
                                                    1,
                                                    0,
                                                    output_file,
                                                    f'_file_{i+1}',
                                                    None)
                out_df = pd.read_csv(output_file, delimiter='\t', dtype=str, index_col=False)
                cols = list(out_df.columns)
                out_df = out_df[[cols[1], cols[0]] + cols[2:]]
                out_df.to_csv(output_file, sep='\t', index=False)
            # Compare the second column of the first file and the first one of the second file
            elif best_match == 'f1s0':
                prf.print_message('Second and first columns are a match', 'info')
                mismatched_f1s0 = old_1_filtered.iloc[:, 1][~old_1_filtered.iloc[:, 1].isin(new_0_filtered.iloc[:, 0])]
                prf.print_message(f"Ids from the consolidated file that didn't match: {mismatched_f1s0}", 'info')
                mismatched_s0f1 = new_0_filtered.iloc[:, 0][~new_0_filtered.iloc[:, 0].isin(old_1_filtered.iloc[:, 1])]
                prf.print_message(f"Ids from the {i+1} file that didn't match: {mismatched_s0f1}", 'info')
                upf.merge_files_by_column_values_df(old_df,
                                                    new_0_filtered,
                                                    1,
                                                    0,
                                                    output_file,
                                                    None,
                                                    f'_file_{i+1}')
                out_df = pd.read_csv(output_file, delimiter='\t', dtype=str, index_col=False)
                cols = list(out_df.columns)
                out_df = out_df[[cols[1], cols[0]] + cols[2:]]
                out_df.to_csv(output_file, sep='\t', index=False)
            # Compare the second columns of both files
            elif best_match == 'f1s1':
                prf.print_message('Second columns are a match', 'info')
                mismatched_f1s1 = old_1_filtered.iloc[:, 1][~old_1_filtered.iloc[:, 1].isin(new_1_filtered.iloc[:, 1])]
                prf.print_message(f"Ids from the consolidated file that didn't match: {mismatched_f1s1}", 'info')
                mismatched_s1f1 = new_1_filtered.iloc[:, 0][~new_1_filtered.iloc[:, 1].isin(old_1_filtered.iloc[:, 1])]
                prf.print_message(f"Ids from the {i+1} file that didn't match: {mismatched_s1f1}", 'info')
                upf.merge_files_by_column_values_df(old_df,
                                                    new_1_filtered,
                                                    1,
                                                    1,
                                                    output_file,
                                                    None,
                                                    f'_file_{i+1}')
                out_df = pd.read_csv(output_file, delimiter='\t', dtype=str, index_col=False)
                cols = list(out_df.columns)
                out_df = out_df[[cols[1], cols[0]] + cols[2:]]
                out_df.to_csv(output_file, sep='\t', index=False)
            else:
                prf.print_message('No columns are a match', 'info')

    # If cleanup TRUE dedeup the final file
    if cleanup:
        prf.print_message('Deduplicating the final file...', 'info')
        ann_df = pd.read_csv(output_file, delimiter='\t', dtype=str, index_col=False)
        ann_df.sort_values(by=ann_df.columns[0], ascending=False, inplace=True)
        ann_df.drop_duplicates(subset=ann_df.columns[0], inplace=True)
        ann_df.to_csv(output_file, sep="\t", index=False)

    return output_file