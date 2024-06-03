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

def get_paths_in_directory(directory):
    """
    Get all paths of files in the specified directory.
    
    Parameters
    ----------
    directory : str
        Path to the directory.
        
    Returns
    -------
    file_paths_dict : list
        List that contains all of the file paths as values.
    """
    # List all files and directories in the specified directory
    all_items = os.listdir(directory)
    
    # Filter out only the paths of files (not directories)
    file_paths = [os.path.join(directory, item) for item in all_items if os.path.isfile(os.path.join(directory, item))]
    
    return file_paths

def get_file_paths_dict(directory):
    """
    Get a dictionary where keys are filenames and values are file paths within the directory.
    
    Parameters
    ----------
    directory : str
        Path to the directory.
        
    Returns
    -------
    file_paths_dict : dict
        Dict that contains all of the filenames as keys and file paths as values.
    """
    file_paths_dict = {}
    
    # List all files and directories in the specified directory
    all_items = os.listdir(directory)
    
    # Iterate over the filenames and construct the dictionary
    for filename in all_items:
        # Construct the file path
        file_path = os.path.join(directory, filename)
        # Check if it's a file (not a directory)
        if os.path.isfile(file_path):
            # Add to the dictionary
            file_paths_dict[filename] = file_path
    
    return file_paths_dict

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