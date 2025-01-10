import os
import shutil
import pandas as pd
from functools import reduce
from typing import Dict, Any, List

def dict_to_df(dictionary: Dict[str, Any]) -> pd.DataFrame:
    """
    Convert a dictionary to a pandas DataFrame.

    Parameters
    ----------
    dictionary : dict
        The dictionary to convert. Keys are used as column headers.

    Returns
    -------
    pd.DataFrame
        The resulting DataFrame, where each key-value pair in the dictionary corresponds to a column in the DataFrame.
    """
    return pd.DataFrame.from_dict(dictionary)


def merge_files_into_same_file_by_key(files: List[str], key_to_merge: str, output_file: str) -> pd.DataFrame:
    """
    Merge multiple TSV files into a single file based on a common key.

    This function reads multiple TSV files, merges them into a single DataFrame based on a common key,
    and writes the merged DataFrame to an output TSV file.

    Parameters
    ----------
    files : List[str]
        List of file paths to the TSV files to be merged.
    key_to_merge : str
        The key column name to merge the files on.
    output_file : str
        The path to the output TSV file where the merged DataFrame will be saved.

    Returns
    -------
    pd.DataFrame
        The merged DataFrame.
    """
    if len(files) == 1:
        shutil.copy(files[0], output_file)
        return None

    # Read all TSV files into a list of DataFrames
    dfs: List[pd.DataFrame] = []
    for file in files:
        current_df: pd.DataFrame = pd.read_csv(file, delimiter='\t', dtype=str, index_col=False)
        dfs.append(current_df)

    # Merge all dataframes based on the key with custom suffixes
    merged_table: pd.DataFrame = dfs[0]
    for i in range(1, len(dfs)):
        suffix = f"_{os.path.basename(files[i]).split('.')[0]}"
        merged_table = pd.merge(merged_table, dfs[i], on=key_to_merge, how='outer', suffixes=('', suffix)).fillna('NA')
    
    # Save the merged table to a TSV file
    merged_table.to_csv(output_file, sep='\t', index=False)

def process_tsv_with_priority(input_file: str, priority_dict: Dict[str, List[str]],
                              output_file: str, best_annotations_bsr: float, 
                              output_columns: List[str]) -> None:
    """
    Process a TSV file based on a dictionary specifying columns and their corresponding values,
    and write the processed DataFrame to a TSV file.

    Parameters
    ----------
    input_file : str
        The path to the input TSV file.
    priority_dict : Dict[str, List[str]]
        A dictionary where the key is a column name (ordered by priority) and the value is a list of column names to return.
    output_file : str
        The path to the output TSV file where the processed DataFrame will be saved.
    best_annotations_bsr : float
        The minimum value threshold for the priority columns.
    output_columns : List[str]
        The list of columns to include in the output TSV file.

    Returns
    -------
    None
    """
    # Read the input TSV file into a DataFrame
    df = pd.read_csv(input_file, delimiter='\t', dtype=str, index_col=False)

    # Initialize a list to store the result rows
    result_rows = []

    # Iterate over each row in the DataFrame
    for index, row in df.iterrows():
        # Initialize a variable to store the selected columns for the current row
        selected_columns = []

        # Iterate over the priority dictionary
        for priority_column, columns_to_return in priority_dict.items():
            # Check if the priority column value is greater than the other priority columns
            if all(float(row[priority_column]) > float(row[other_column]) for other_column in priority_dict if other_column != priority_column):
                # Check if the minimum value in the priority columns meets the best_annotations_bsr threshold
                if float(row[priority_column]) >= best_annotations_bsr:
                    selected_columns.extend(columns_to_return)
                    break

        # If selected_columns is not empty, extract the values for the selected columns
        if selected_columns:
            # Get database name
            db_name = selected_columns[0].split('_')[0]
            # Get the selected columns from the row
            result_row = row[selected_columns]
        else:
            # If the value is smaller than the best_annotations_bsr, get the 'Locus' column
            result_row = row[['Locus']]

        # Append the result row to the list
        result_rows.append(result_row)

    # Concatenate all result rows into a single DataFrame
    result_df = pd.concat(result_rows, axis=1).T

    # Write the processed DataFrame to a TSV file with specific columns
    result_df.to_csv(output_file, sep='\t', index=False, columns=output_columns)
