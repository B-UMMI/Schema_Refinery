import os
import shutil
import pandas as pd
from typing import Dict, Any, List, Union

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

def merge_files_by_column_values(file1: str, file2: str, column_value1: Union[str, int], column_value2: Union[str, int], output_file: str) -> pd.DataFrame:
    """
    Merge two TSV files into a single file based on specified column values or indices.

    This function reads two TSV files, merges them into a single DataFrame based on specified column values or indices,
    and writes the merged DataFrame to an output TSV file.

    Parameters
    ----------
    file1 : str
        File path to the first TSV file.
    file2 : str
        File path to the second TSV file.
    column_value1 : Union[str, int]
        The column value or index to merge the first file on.
    column_value2 : Union[str, int]
        The column value or index to merge the second file on.
    output_file : str
        The path to the output TSV file where the merged DataFrame will be saved.

    Returns
    -------
    pd.DataFrame
        The merged DataFrame.
    """
    # Read the TSV files into DataFrames
    df1 = pd.read_csv(file1, delimiter='\t', dtype=str, index_col=False)
    df2 = pd.read_csv(file2, delimiter='\t', dtype=str, index_col=False)

    # Convert column indices to column names if necessary
    if isinstance(column_value1, int):
        column_value1 = df1.columns[column_value1]
    if isinstance(column_value2, int):
        column_value2 = df2.columns[column_value2]

    # Merge the dataframes based on the specified column values
    merged_table = pd.merge(df1, df2, left_on=column_value1, right_on=column_value2, how='outer')
    
    # Drop the 'Locus_y' which is original subject id column and rename 'Locus_x' to 'Locus'
    if 'Locus_y' in merged_table.columns:
        merged_table.drop(columns=['Locus_y'], inplace=True)
    if 'Locus_x' in merged_table.columns:
        merged_table.rename(columns={'Locus_x': 'Locus'}, inplace=True)

    # Rename specified columns by adding 'matched_' prefix
    columns_to_rename = {
        'Protein_ID': 'matched_Protein_ID',
        'Protein_product': 'matched_Protein_product',
        'Protein_short_name': 'matched_Protein_short_name',
        'Protein_BSR': 'matched_Protein_BSR'
    }
    merged_table.rename(columns=columns_to_rename, inplace=True)

    # Save the merged table to a TSV file
    merged_table.to_csv(output_file, sep='\t', index=False)

    return merged_table

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
            db_name = selected_columns[1].split('_')[0]
            # Get the selected columns from the row as a dictionary
            selected_columns_dict = row[selected_columns].to_dict()
            # Add the database name to the dictionary
            selected_columns_dict['Database'] = db_name

            # Rename the dictionary keys using new_column_names
            selected_columns_dict = {output_columns[i]: v for i, v in enumerate(selected_columns_dict.values())}
            result_row = selected_columns_dict
        else:
            # If the value is smaller than the best_annotations_bsr, get the 'Locus' column
            result_row = row[['Locus']].to_dict()
            # Add a placeholder for the missing columns
            for col in output_columns:
                if col not in result_row:
                    result_row[col] = None

        # Append the result row to the list
        result_rows.append(result_row)

    # Concatenate all result rows into a single DataFrame
    result_df = pd.DataFrame(result_rows)

    # Write the processed DataFrame to a TSV file with specific columns
    result_df.to_csv(output_file, sep='\t', index=False, columns=output_columns)
