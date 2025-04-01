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
		merged_table = pd.merge(merged_table, dfs[i], on=key_to_merge, how='left', suffixes=('', suffix)).fillna('NA')
	
	# Save the merged table to a TSV file
	merged_table.to_csv(output_file, sep='\t', index=False)

def merge_files_by_column_values(file1: str, file2: str, column_value1: Union[str, int], column_value2: Union[str, int], output_file: str) -> str:
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
	str
		The path to the tsv file of the merged DataFrame.
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
	merged_table = pd.merge(df1, df2, left_on=column_value1, right_on=column_value2, how='left')
	
	# Drop all the locus columns
	if 'Locus_y' in merged_table.columns:
		merged_table.drop(columns=['Locus_y'], inplace=True)
	if 'Locus_x' in merged_table.columns:
		merged_table.drop(columns=['Locus_x'], inplace=True)
	if 'Query' in merged_table.columns:
		if 'Locus' in merged_table.columns:
			merged_table.drop(columns=['Locus'], inplace=True)

	# Rename specified columns by adding 'matched_' prefix
	columns_to_rename = {
		'Proteome_ID': 'matched_Proteome_ID',
		'Proteome_product': 'matched_Proteome_product',
		'Proteome_gene_name': 'matched_Proteome_gene_name',
		'Proteome_BSR': 'matched_Proteome_BSR',
		'Proteome_ID_best_proteomes_annotations_swiss_prot': 'matched_Proteome_ID_best_proteomes_annotations_swiss_prot',
		'Proteome_product_best_proteomes_annotations_swiss_prot': 'matched_Proteome_product_best_proteomes_annotations_swiss_prot',
		'Proteome_gene_name_best_proteomes_annotations_swiss_prot': 'matched_Proteome_gene_name_best_proteomes_annotations_swiss_prot',
		'Proteome_BSR_best_proteomes_annotations_swiss_prot': 'matched_Proteome_BSR_best_proteomes_annotations_swiss_prot'

	}
	merged_table.rename(columns=columns_to_rename, inplace=True)

	# Save the merged table to a TSV file
	merged_table.to_csv(output_file, sep='\t', index=False)

	return output_file

def merge_files_by_column_values_df(df1: pd.DataFrame, df2: pd.DataFrame, column_value1: Union[str, int], column_value2: Union[str, int], output_file: str, left: str, right: str) -> pd.DataFrame:
	"""
	Merge two TSV files into a single file based on specified column values or indices.

	This function reads two pandas DataFrames, merges them into a single DataFrame based on specified column values or indices,
	and writes the merged DataFrame to an output TSV file.

	Parameters
	----------
	df1 : pd.DataFrame
		First DataFrame.
	df2 : pd.DataFrame
		Second DataFrame
	column_value1 : Union[str, int]
		The column value or index to merge the first file on.
	column_value2 : Union[str, int]
		The column value or index to merge the second file on.
	output_file : str
		The path to the output TSV file where the merged DataFrame will be saved.
	left : str
		sufix of the left DataFrame
	right : str
		sufix of the right DataFrame

	Returns
	-------
	pd.DataFrame
		The merged DataFrame.
	"""
	# Convert column indices to column names if necessary
	if isinstance(column_value1, int):
		column_value1 = df1.columns[column_value1]
	if isinstance(column_value2, int):
		column_value2 = df2.columns[column_value2]

	# Merge the dataframes based on the specified column values
	merged_table = pd.merge(df1, df2, left_on=column_value1, right_on=column_value2, how='left', suffixes=(left, right))
	
	# Drop the 'Locus_y' which is original subject id column and rename 'Locus_x' to 'Locus'
	# Drop all the locus columns
	if 'Locus_y' in merged_table.columns:
		merged_table.drop(columns=['Locus_y'], inplace=True)
	if 'Locus_x' in merged_table.columns:
		merged_table.drop(columns=['Locus_x'], inplace=True)
	if 'Query' in merged_table.columns:
		if 'Locus' in merged_table.columns:
			merged_table.drop(columns=['Locus'], inplace=True)

	# Save the merged table to a TSV file
	merged_table.to_csv(output_file, sep='\t', index=False)

	return merged_table
