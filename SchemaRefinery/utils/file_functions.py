import os
import shutil
import pandas as pd
import csv
from itertools import zip_longest, islice
from typing import List, Dict, Any, Union

try:
    from utils import (iterable_functions as itf)
except ModuleNotFoundError:
    from SchemaRefinery.utils import (iterable_functions as itf)

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

def check_and_delete_file(file: str) -> None:
    """
    Deletes a file based on the input file path.

    Parameters
    ----------
    file : str
        File path.

    Returns
    -------
    None
    """
    if os.path.isfile(file):
        os.remove(file)

def copy_file(source_file: str, destination_file: str) -> None:
    """
    Copies the source file to the destination.

    Parameters
    ----------
    source_file : str
        Path to the source file.
    destination_file : str
        Path to where to copy the file.

    Returns
    -------
    None
    """
    shutil.copy(source_file, destination_file)

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

def get_paths_dict(directory: str, type_: str) -> Dict[str, str]:
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
    Dict[str, str]
        A dictionary with filenames as keys and their full paths as values. The contents
        are filtered based on the `type_` parameter.

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

    paths_dict: Dict[str, str] = {}
    all_items: List[str] = os.listdir(directory)
    
    for filename in all_items:
        file_path: str = os.path.join(directory, filename)
        if if_type(file_path):
            paths_dict[filename] = file_path
    
    return paths_dict

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

def write_dict_to_tsv(file_path: str, data: Dict[str, List[Any]]) -> None:
    """
    Write a dictionary to a TSV (Tab-Separated Values) file.

    Parameters
    ----------
    file_path : str
        The file path to save the TSV file.
    data : Dict[str, List[Any]]
        The dictionary where keys are column names and values are lists of values for each column.

    Returns
    -------
    None
    """
    with open(file_path, 'w') as f:
        f.write('\t'.join(data.keys()) + '\n')
        for row in zip_longest(*data.values(), fillvalue=''):
            row_str: str = '\t'.join(map(str, row))
            f.write(row_str + '\n')

def concat_files(source_file: str, destination_file: str) -> None:
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

def merge_folders(folder1: str, folder2: str, output_folder: str) -> None:
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

    Returns
    -------
    None

    Notes
    -----
    - If the output folder does not exist, it will be created.
    - If there are file naming conflicts, the conflicting files will be renamed with a "_copy" suffix.
    - The function ensures that all intermediate directories are created as needed.
    """
    def copy_files(src_folder: str) -> None:
        """
        Copy files from the source folder to the output folder.

        Parameters
        ----------
        src_folder : str
            Path to the source folder from which files will be copied.

        Returns
        -------
        None
        """
        for root, _, files in os.walk(src_folder):
            for file in files:
                src_file_path: str = os.path.join(root, file)
                relative_path: str = os.path.relpath(src_file_path, src_folder)
                dest_file_path: str = os.path.join(output_folder, relative_path)
                
                os.makedirs(os.path.dirname(dest_file_path), exist_ok=True)
                
                if os.path.exists(dest_file_path):
                    base, ext = os.path.splitext(dest_file_path)
                    dest_file_path = f"{base}_copy{ext}"
                
                shutil.copy2(src_file_path, dest_file_path)
    
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)
    
    copy_files(folder1)
    copy_files(folder2)