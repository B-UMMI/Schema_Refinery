import os
import sys
import shutil
import pandas as pd
import csv
from itertools import islice
from typing import List, Dict, Any, Union, Tuple

try:
	from utils import (sequence_functions as sf,
					   print_functions as pf)
	from AdaptLoci import AdaptLoci as al
except ModuleNotFoundError:
	from SchemaRefinery.utils import (sequence_functions as sf,
									  print_functions as pf)
	from SchemaRefinery.AdaptLoci import AdaptLoci as al

def create_directory(dir: str) -> bool:
	"""
	Creates a directory based on the input dir path.

	Parameters
	----------
	dir : str
		Directory path.

	Returns
	-------
	bool
		True if the directory was created, False if it already exists.
	"""
	if not os.path.exists(dir):
		os.makedirs(dir)
		return True
	else:
		return False


def import_df_from_file(file_path: str, sep: str) -> pd.DataFrame:
	"""
	Using pandas, imports a file path as a dataframe.

	Parameters
	----------
	file_path : str
		Path to the file.
	sep : str
		String to separate entries in the file.

	Returns
	-------
	pd.DataFrame
		Pandas dataframe that contains the file values separated by sep.
	"""
	df: pd.DataFrame = pd.read_csv(file_path, sep=sep)
	return df


def get_paths_in_directory(directory: str, type_: str) -> List[str]:
	"""
	Get all paths of items in the specified directory filtered by type.

	This function lists all items (files, directories, or both) in the given directory
	and returns their full paths. The items to be listed can be filtered by specifying
	a type: 'files' for files only, 'directories' for directories only, or 'all' for both.

	Parameters
	----------
	directory : str
		The path to the directory from which to retrieve item paths.
	type_ : str
		The type of items to include in the output list. Valid options are:
		- 'files': Include only files.
		- 'directories': Include only directories.
		- 'all': Include both files and directories.

	Returns
	-------
	List[str]
		A list containing the full paths of the items in the directory, filtered by the specified type.

	Raises
	------
	ValueError
		If the `type_` parameter is not one of the valid options ('files', 'directories', 'all').
	"""
	if type_ == 'files':
		if_type = os.path.isfile
	elif type_ == 'directories':
		if_type = os.path.isdir
	elif type_ == 'all':
		if_type = os.path.exists
	else:
		raise ValueError(f"Invalid type: {type_}")

	all_items: List[str] = os.listdir(directory)
	file_paths: List[str] = [os.path.join(directory, item) for item in all_items if if_type(os.path.join(directory, item))]
	
	return file_paths


def get_paths_in_directory_with_suffix(directory: str, suffix: str) -> List[str]:
	"""
	Get all paths of files in the specified directory that end with a given suffix.

	Parameters
	----------
	directory : str
		Path to the directory.
	suffix : str
		The suffix that the files must end with.

	Returns
	-------
	List[str]
		List that contains all of the file paths with the specified suffix.
	"""
	all_items: List[str] = os.listdir(directory)
	file_paths: List[str] = [os.path.join(directory, item) for item in all_items if os.path.isfile(os.path.join(directory, item)) and item.endswith(suffix)]
	
	return file_paths


def map_basename_to_path(paths, extension=False, remove_suffix=None) -> List[str]:
	"""
	Map file basenames to their full paths.

	Parameters
	----------
	paths : list
		List with the full paths to the files.

	Returns
	-------
	mapped_basenames : dict
		A dictionary with the file basenames and their full paths.
			- Keys are the basenames of the files.
			- Values are the full paths to the files.
	"""
	mapped_basenames: Dict[str, str] = {}
	for path in paths:
		basename: str = os.path.basename(path)
		if not extension:
			basename = '.'.join(basename.split('.')[:-1])
		if remove_suffix:
			basename = basename.replace(remove_suffix, '')
		mapped_basenames[basename] = path

	return mapped_basenames


def file_basename(file_path: str, file_extension: bool = True) -> str:
	"""
	Get the file name from a file path.

	Parameters
	----------
	file_path : str
		The path to the file.
	file_extension : bool, optional
		Whether to include the file extension in the returned name (default is True).

	Returns
	-------
	str
		The file name extracted from the file path.
	"""
	if not file_extension:
		return os.path.basename(file_path).split('.')[0]
	else:
		return os.path.basename(file_path)


def join_paths(parent_path: str, child_paths: List[str]) -> str:
	"""
	Create a path by joining a parent directory and a list of child paths.

	Parameters
	----------
	parent_path : str
		The parent directory path.
	child_paths : List[str]
		A list of child paths to join with the parent path.

	Returns
	-------
	str
		The joined path.
	"""
	joined_paths: str = os.path.join(parent_path, *child_paths)
	return joined_paths


def read_lines(input_file: str, strip: bool = True, num_lines: Union[int, None] = None) -> List[str]:
	"""
	Read lines in a file.

	Parameters
	----------
	input_file : str
		Path to the input file.
	strip : bool, optional
		Specify if lines should be stripped of leading and trailing white spaces and new line characters (default is True).
	num_lines : Union[int, None], optional
		Number of lines to read from the file (default is None, which reads all lines).

	Returns
	-------
	List[str]
		List with the lines read from the input file.
	"""
	with open(input_file, 'r') as infile:
		if num_lines is None:
			lines: List[str] = [line for line in infile.readlines()]
		else:
			lines = list(islice(infile, num_lines))

	if strip:
		lines = [line.strip() for line in lines]

	return lines


def write_lines(lines: List[str], output_file: str, joiner: str = '\n', write_mode: str = 'w') -> None:
	"""
	Write a list of strings to a file.

	Parameters
	----------
	lines : List[str]
		List with the lines/strings to write to the output file.
	output_file : str
		Path to the output file.
	joiner : str, optional
		Character used to join lines (default is '\n').
	write_mode : str, optional
		Specify write mode ('w' creates file if it does not exist and truncates and over-writes existing file,
		'a' creates file if it does not exist and appends to the end of file if it exists) (default is 'w').

	Returns
	-------
	None
	"""
	joined_lines: str = joiner.join(lines)
	write_to_file(joined_lines, output_file, write_mode, '\n')


def write_to_file(text: str, output_file: str, write_mode: str, end_char: str) -> None:
	"""
	Write a single string to a file.

	Parameters
	----------
	text : str
		A single string to write to the output file.
	output_file : str
		Path to the output file.
	write_mode : str
		Specify write mode ('w' creates file if it does not exist and truncates and over-writes existing file,
		'a' creates file if it does not exist and appends to the end of file if it exists).
	end_char : str
		Character added to the end of the file.

	Returns
	-------
	None
	"""
	with open(output_file, write_mode) as out:
		out.write(text + end_char)


def read_tabular(input_file: str, delimiter: str = '\t') -> List[List[str]]:
	"""
	Read a tabular (TSV) file.

	Parameters
	----------
	input_file : str
		Path to a tabular file.
	delimiter : str, optional
		Delimiter used to separate file fields (default is '\t').

	Returns
	-------
	List[List[str]]
		A list with a sublist per line in the input file. Each sublist has the fields that were separated by the defined delimiter.
	"""
	with open(input_file, 'r') as infile:
		reader = csv.reader(infile, delimiter=delimiter)
		lines: List[List[str]] = [line for line in reader]

	return lines


def copy_folder(src_folder: str, dest_folder: str) -> None:
	"""
	Copy the contents of one folder to another folder.

	Parameters
	----------
	src_folder : str
		Path to the source folder to copy.
	dest_folder : str
		Path to the destination folder where the contents will be copied.

	Returns
	-------
	None

	Notes
	-----
	- If the destination folder does not exist, it will be created.
	- The function copies the entire folder, including all subdirectories and files.
	- If the destination folder already exists, existing files will be overwritten.
	"""
	if not os.path.exists(dest_folder):
		os.makedirs(dest_folder)
	
	shutil.copytree(src_folder, dest_folder, dirs_exist_ok=True)


def merge_folders(folder1: str, folder2: str, output_folder: str, cpu: int, bsr: float, translation_table: int) -> None:
	"""
	Merge the contents of two folders into an output folder.

	Parameters
	----------
	folder1 : str
		Path to the first folder to merge.
	folder2 : str
		Path to the second folder to merge.
	output_folder : str
		Path to the output folder where the merged contents will be stored.
	cpu : int
        Number of CPU cores that will be used to run the process.
    blast_score_ratio : float
        The BLAST Score Ratio value that will be used to evaluate
        allele similarity and select representative alleles.
    translation_table : int
        Genetic code used to translate alleles.

	Returns
	-------
	None

	Notes
	-----
	- The function ensures that all intermediate directories are created as needed.
	"""
	def get_fasta_files(folder: str) ->  List[str]:
		"""
		Get a list of .fasta files in the given folder.

		Parameters
		----------
		folder : str
			The path to the folder to search for .fasta files.

		Returns
		-------
		List[str]
			A list of .fasta file names in the folder.
		"""
		folder_fastas: List[str] = []

		folder_fastas = [f for f in os.listdir(folder) if os.path.isfile(os.path.join(folder, f)) and f.endswith('.fasta')]

		return folder_fastas
	def read_and_merge_fasta_files(folder1: str, folder2: str, locus: str) -> Dict[int, str]:
		"""
		Read and merge .fasta files from two folders for a given locus.

		Parameters
		----------
		folder1 : str
			Path to the first folder.
		folder2 : str
			Path to the second folder.
		locus : str
			The name of the locus to merge.

		Returns
		-------
		Dict[int, str]
			A dictionary with merged .fasta sequences, with keys as sequential integers starting from 1.
		"""
		folder1_locus = os.path.join(folder1, locus)
		folder2_locus = os.path.join(folder2, locus)

		folder1_locus_fasta = sf.read_fasta_file_dict(folder1_locus)
		folder2_locus_fasta = sf.read_fasta_file_dict(folder2_locus)

		merged_locus_fasta = {**folder1_locus_fasta, **folder2_locus_fasta}
		
		seen_hashes = set()
		merged_locus_fasta = {i: str(value.seq) for i, (key, value) in enumerate(
			{k: v for k, v in merged_locus_fasta.items() if sf.hash_sequence(str(v)) not in seen_hashes and not seen_hashes.add(sf.hash_sequence(str(v)))}.items(), 1)}
		return merged_locus_fasta

	def write_fasta_file(fasta_dict: Dict[int, str], output_path: str) -> None:
		"""
		Write a dictionary of .fasta sequences to a file.

		Parameters
		----------
		fasta_dict : Dict[int, str]
			A dictionary with .fasta sequences.
		output_path : str
			The path to the output .fasta file.
		"""
		with open(output_path, 'w') as fasta_file:
			for allele, sequence in fasta_dict.items():
				fasta_file.write(f">{allele}\n{sequence}\n")
			fasta_file.write("\n")

	def process_loci_with_same_name(folder1: str, folder2: str, loci_with_same_name: set, temp_folder: str) -> None:
		"""
		Process loci with the same name from two folders and merge them.

		Parameters
		----------
		folder1 : str
			Path to the first folder.
		folder2 : str
			Path to the second folder.
		loci_with_same_name : set
			A set of loci names that are present in both folders.
		temp_folder : str
			Path to the temporary folder to store merged loci.
		"""
		for locus in loci_with_same_name:
			merged_locus_fasta = read_and_merge_fasta_files(folder1, folder2, locus)
			output_locus = os.path.join(temp_folder, locus)
			write_fasta_file(merged_locus_fasta, output_locus)

	def process_remaining_files(files: List[str], folder: str, temp_folder: str) -> None:
		"""
		Process remaining .fasta files from a folder and copy them to a temporary folder.

		Parameters
		----------
		files : List[str]
			A list of .fasta file names.
		folder : str
			Path to the folder containing the .fasta files.
		temp_folder : str
			Path to the temporary folder to store the processed files.
		"""
		for file in files:
			temp_fasta = os.path.join(temp_folder, os.path.basename(file))
			fasta_dict = sf.read_fasta_file_dict(file)
			fasta_dict = {i: str(value.seq) for i, (key, value) in enumerate(fasta_dict.items(), 1)}
			write_fasta_file(fasta_dict, temp_fasta)

	# Create a temporary folder to store intermediate files
	temp_folder = os.path.join(output_folder, 'temp')
	create_directory(temp_folder)

	# Get the list of .fasta files from both folders
	folder1_files = get_fasta_files(folder1)
	folder2_files = get_fasta_files(folder2)

	# Find loci with the same name in both folders
	loci_with_same_name = set(folder1_files) & set(folder2_files)
	if loci_with_same_name:
		pf.print_message(f"The following loci have the same name in both folders and will be merged: {loci_with_same_name}", "warning")
		process_loci_with_same_name(folder1, folder2, loci_with_same_name, temp_folder)
	else:
		pf.print_message(f"There is no overlap between the sets of loci.", "warning")

	# Remove loci with the same name from the lists
	folder1_files = [f for f in folder1_files if f not in loci_with_same_name]
	folder2_files = [f for f in folder2_files if f not in loci_with_same_name]

	# Add full paths to the remaining files
	folder1_files = [os.path.join(folder1, f) for f in folder1_files]
	folder2_files = [os.path.join(folder2, f) for f in folder2_files]

	# Process the remaining files and copy them to the temporary folder
	process_remaining_files(folder1_files, folder1, temp_folder)
	process_remaining_files(folder2_files, folder2, temp_folder)

	# Write the AdaptLoci.txt file and run the AdaptLoci process
	pf.print_message("Adapting loci from both folders into one Schema")
	temp_output = os.path.join(output_folder, 'temp_out')
	# Run the adapt loci module
	al.adapt_loci(temp_folder, temp_output, cpu, bsr, translation_table)

	# Fix IDs in the .fasta files in the output folder
	copy_folder(temp_output, output_folder)

	# Remove the temporary folder
	shutil.rmtree(temp_folder)
	shutil.rmtree(temp_output)

def cleanup(directory: str, exclude: List[str]) -> None:
	"""
	Clean up a directory by removing all files and subdirectories except those specified in the exclusion list.

	Parameters
	----------
	directory : str
		The path to the directory to clean up.
	exclude : List[str]
		A list of filenames or subdirectory names to exclude from removal.

	Returns
	-------
	None
	"""
	for item in os.listdir(directory):
		item_path = os.path.join(directory, item)
		if item_path not in exclude:
			if os.path.isfile(item_path) or os.path.islink(item_path):
				os.remove(item_path)
			elif os.path.isdir(item_path):
				shutil.rmtree(item_path)


def concatenate_files(files, output_file, header=None):
	"""Concatenate files.

	Parameters
	----------
	files : list
		List with the paths to the files to concatenate.
	output_file : str
		Path to the output file that will store the concatenation
		of input files.
	header : str or NoneType
		Specify a header that should be written as the first line
		in the output file.

	Returns
	-------
	output_file : str
		Path to the output file that was created with
		the concatenation of input files.
	"""
	with open(output_file, 'w') as outfile:
		if header is not None:
			outfile.write(header)
		for file in files:
			with open(file, 'r') as infile:
				shutil.copyfileobj(infile, outfile)

	return output_file
