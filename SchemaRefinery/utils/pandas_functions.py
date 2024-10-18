import os
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

def merge_files_into_same_file_by_key(files: List[str], key_to_merge, output_file) -> pd.DataFrame:
    dfs: List[pd.DataFrame] = []
    for file in files:
        current_df: pd.DataFrame = pd.read_csv(file, delimiter='\t', dtype=str)
        dfs.append(current_df)
    
    # Merge all dataframes based on the key with custom suffixes
    merged_table: pd.DataFrame = dfs[0]
    for i in range(1, len(dfs)):
        suffix = f"_{os.path.basename(files[i]).split('.')[0]}"
        merged_table = pd.merge(merged_table, dfs[i], on=[key_to_merge], how='left', suffixes=('', suffix)).fillna('NA')
    

    # Save the merged table to a TSV file
    merged_table.to_csv(output_file, sep='\t', index=False)

    return output_file