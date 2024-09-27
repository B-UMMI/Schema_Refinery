import os
import shutil
import pandas as pd
import csv
from itertools import zip_longest, islice

try:
    from utils import (iterable_functions as itf)
except ModuleNotFoundError:
    from SchemaRefinery.utils import (iterable_functions as itf)

def create_directory(dir):
    """
    Creates directory based on input dir path.

    Parameters
    ----------
    file : str
        dir path

    Returns
    -------
    return : None
        Creates directory at desired dir path.
    """

    if not os.path.exists(dir):
        os.makedirs(dir)
        return True
    else:
        return False

def check_and_delete_file(file:str):
    """
    Deletes file based on input file path.

    Parameters
    ----------
    file : str
        File path

    Returns
    -------
    return : None
        Deletes file based on input file path
    """
    
    if os.path.isfile(file):
        os.remove(file)

def copy_file(source_file, destination_file):
    """
    Copies source file to the destination.
    
    Parameters
    ----------
    source_file : str
        Path to the file.
    destination_file : str
        Path to where to copy the file.
    
    Returns
    -------
    return : None
        Copies the file.
    """
    # Copy the file to the destination directory
    shutil.copy(source_file, destination_file)

def import_df_from_file(file_path, sep):
    """
    Using pandas imports an file path as dataframe.
    
    Parameters
    ----------
    file_path : str
        Path to the file.
    sep : str
        By which string to seprate entries in the file.

    Returns
    -------
    df : pandas dataframe
        Pandas dataframe that contains the file values seperated by sep.
    """
    df = pd.read_csv(file_path, sep=sep)

    return df

def get_paths_in_directory(directory, type_):
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
    list
        A list containing the full paths of the items in the directory, filtered by the specified type.

    Raises
    ------
    ValueError
        If the `type_` parameter is not one of the valid options ('files', 'directories', 'all').

    Examples
    --------
    >>> get_paths_in_directory('/path/to/directory', 'files')
    ['/path/to/directory/file1.txt', '/path/to/directory/file2.jpg']

    >>> get_paths_in_directory('/path/to/directory', 'directories')
    ['/path/to/directory/subdir1', '/path/to/directory/subdir2']

    >>> get_paths_in_directory('/path/to/directory', 'all')
    ['/path/to/directory/file1.txt', '/path/to/directory/subdir1']
    """
    if type_ == 'files':
        if_type = os.path.isfile
    elif type_ == 'directories':
        if_type = os.path.isdir
    elif type_ == 'all':
        if_type = os.path.exists
    else:
        raise ValueError(f"Invalid type: {type_}")
    # List all files and directories in the specified directory
    all_items = os.listdir(directory)
    
    # Filter out only the paths of files (not directories)
    file_paths = [os.path.join(directory, item) for item in all_items if if_type(os.path.join(directory, item))]
    
    return file_paths

def get_paths_dict(directory, type_):
    """
    Get a dictionary where keys are filenames and values are file paths within the directory,
    filtered by the specified type: files, directories, or all items.

    Parameters
    ----------
    directory : str
        The path to the directory from which to retrieve file paths.
    type_ : str
        The type of items to include in the output dictionary. Valid options are:
        - 'files': Include only files.
        - 'directories': Include only directories.
        - 'all': Include both files and directories.

    Returns
    -------
    dict
        A dictionary with filenames as keys and their full paths as values. The contents
        are filtered based on the `type_` parameter.

    Raises
    ------
    ValueError
        If the `type_` parameter is not one of the valid options ('files', 'directories', 'all').

    Examples
    --------
    >>> get_paths_dict('/path/to/directory', 'files')
    {'file1.txt': '/path/to/directory/file1.txt', 'file2.txt': '/path/to/directory/file2.txt'}

    >>> get_paths_dict('/path/to/directory', 'directories')
    {'subdir1': '/path/to/directory/subdir1', 'subdir2': '/path/to/directory/subdir2'}

    >>> get_paths_dict('/path/to/directory', 'all')
    {'file1.txt': '/path/to/directory/file1.txt', 'subdir1': '/path/to/directory/subdir1'}
    """

    if type_ == 'files':
        if_type = os.path.isfile
    elif type_ == 'directories':
        if_type = os.path.isdir
    elif type_ == 'all':
        if_type = os.path.exists
    else:
        raise ValueError(f"Invalid type: {type_}")
    paths_dict = {}
    
    # List all files and directories in the specified directory
    all_items = os.listdir(directory)
    
    # Iterate over the filenames and construct the dictionary
    for filename in all_items:
        # Construct the file path
        file_path = os.path.join(directory, filename)
        # Check if it's a file (not a directory)
        if if_type(file_path):
            # Add to the dictionary
            paths_dict[filename] = file_path
    
    return paths_dict

def get_paths_in_directory_with_suffix(directory, suffix):
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
    file_paths : list
        List that contains all of the file paths with the specified suffix.
    """
    # List all files and directories in the specified directory
    all_items = os.listdir(directory)
    
    # Filter out only the paths of files (not directories) that end with the specified suffix
    file_paths = [os.path.join(directory, item) for item in all_items if os.path.isfile(os.path.join(directory, item)) and item.endswith(suffix)]
    
    return file_paths

def write_dict_to_tsv(file_path, data):
    """
    Write a dictionary to a TSV (Tab-Separated Values) file.
    
    Parameters
    ----------
    file_path : str
        The file path to save the TSV file.
    data : dict 
        The dictionary where keys are column names and values are lists of values for each column.
    
    Returns
    -------
    return : None
        Write dict to file.
    """
    with open(file_path, 'w') as f:
        # Write headers
        f.write('\t'.join(data.keys()) + '\n')

        # Write data
        for row in zip_longest(*data.values(), fillvalue=''):
            row_str = '\t'.join(map(str, row))
            f.write(row_str + '\n')

def concat_files(source_file, destination_file):
    """
    Concatenates the source file to the destination file.
    
    Parameters
    ----------
    source_file : str
        The path to the source file.
    destination_file : str
        The path to the destination file.

    Returns
    -------
    None
    """
    with open(destination_file, 'a') as outfile, open(source_file, 'r') as infile:
        shutil.copyfileobj(infile, outfile)
        
def file_basename(file_path, file_extension=True):
    """
    Get the file name from a file path.
    
    Parameters
    ----------
    file_path : str
        The path to the file.
    
    Returns
    -------
    str
        The file name extracted from the file path.
    """
    if not file_extension:
        return os.path.basename(file_path).split('.')[0]
    else:
        return os.path.basename(file_path)

def join_paths(parent_path, child_paths):
	"""Create path by joining a parent directory and a list of child paths."""
	joined_paths = os.path.join(parent_path, *child_paths)

	return joined_paths

def read_lines(input_file, strip=True, num_lines=None):
	"""Read lines in a file.

	Parameters
	----------
	input_file : str
		Path to the input file.
	strip : bool
		Specify if lines should be stripped of leading
		and trailing white spaces and new line characters.

	Returns
	-------
	lines : list
		List with the lines read from the input file.
	"""
	with open(input_file, 'r') as infile:
		if num_lines is None:
			lines = [line for line in infile.readlines()]
		else:
			lines = list(islice(infile, num_lines))

	if strip is True:
		lines = [line.strip() for line in lines]

	return lines


def write_lines(lines, output_file, joiner='\n', write_mode='w'):
	"""Write a list of strings to a file.

	Parameters
	----------
	lines : list
		List with the lines/strings to write to the output
		file.
	output_file : str
		Path to the output file.
	joiner : str
		Character used to join lines.
	write_mode : str
		Specify write mode ('w' creates file if it does not
		exist and truncates and over-writes existing file,
		'a' creates file if it does not exist and appends to
		the end of file if it exists).
	"""
	joined_lines = itf.join_list(lines, joiner)

	write_to_file(joined_lines, output_file, write_mode, '\n')

def write_to_file(text, output_file, write_mode, end_char):
	"""Write a single string to a file.

	Parameters
	----------
	text : str
		A single string to write to the output file.
	output_file : str
		Path to the output file.
	write_mode : str
		Specify write mode ('w' creates file if it does not
		exist and truncates and over-writes existing file,
		'a' creates file if it does not exist and appends to
		the end of file if it exists.).
	end_char : str
		Character added to the end of the file.
	"""
	with open(output_file, write_mode) as out:
		out.write(text+end_char)
  
def read_tabular(input_file, delimiter='\t'):
	"""Read a tabular (TSV) file.

	Parameters
	----------
	input_file : str
		Path to a tabular file.
	delimiter : str
		Delimiter used to separate file fields.

	Returns
	-------
	lines : list
		A list with a sublist per line in the input file.
		Each sublist has the fields that were separated by
		the defined delimiter.
	"""
	with open(input_file, 'r') as infile:
		reader = csv.reader(infile, delimiter=delimiter)
		lines = [line for line in reader]

	return lines

def copy_folder(src_folder, dest_folder):
    """
    Copy the contents of one folder to another folder.

    Parameters
    ----------
    src_folder : str
        Path to the source folder to copy.
    dest_folder : str
        Path to the destination folder where the contents will be copied.

    Notes
    -----
    - If the destination folder does not exist, it will be created.
    - The function copies the entire folder, including all subdirectories and files.
    - If the destination folder already exists, existing files will be overwritten.
    """
    # Ensure the destination folder exists
    if not os.path.exists(dest_folder):
        os.makedirs(dest_folder)
    
    # Copy the entire folder
    shutil.copytree(src_folder, dest_folder, dirs_exist_ok=True)

def merge_folders(folder1, folder2, output_folder):
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

    Notes
    -----
    - If the output folder does not exist, it will be created.
    - If there are file naming conflicts, the conflicting files will be renamed with a "_copy" suffix.
    - The function ensures that all intermediate directories are created as needed.
    """
    # Helper function to copy files from a folder to the output folder
    def copy_files(src_folder):
        """
        Copy files from the source folder to the output folder.

        Parameters
        ----------
        src_folder : str
            Path to the source folder from which files will be copied.
        """
        for root, _, files in os.walk(src_folder):
            for file in files:
                src_file_path = os.path.join(root, file)
                relative_path = os.path.relpath(src_file_path, src_folder)
                dest_file_path = os.path.join(output_folder, relative_path)
                
                # Ensure the destination directory exists
                os.makedirs(os.path.dirname(dest_file_path), exist_ok=True)
                
                # Handle file naming conflicts
                if os.path.exists(dest_file_path):
                    base, ext = os.path.splitext(dest_file_path)
                    dest_file_path = f"{base}_copy{ext}"
                
                # Copy the file
                shutil.copy2(src_file_path, dest_file_path)
    
        # Ensure the output folder exists
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)
    
    # Copy files from both folders
    copy_files(folder1)
    copy_files(folder2)