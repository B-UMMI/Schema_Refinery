import os
import shutil
import pandas as pd
from itertools import zip_longest

def create_directory(dir:str):
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

    if not os.path.isdir(dir):
        os.mkdir(dir)

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