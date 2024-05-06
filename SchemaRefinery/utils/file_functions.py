import os
import pandas as pd

def create_directory(dir:str):
    """
    Creates directory based on input dir path.

    Parameters
    ----------
    file : str
        dir path

    Returns
    -------
    return : No return
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
    return : No return
        Deletes file based on input file path
    """
    
    if os.path.isfile(file):
        os.remove(file)

def import_df_from_file(file_path, sep):
    df = pd.read_csv(file_path, sep=sep)

    return df

def get_paths_in_directory(directory):
    """
    Get all paths of files in the specified directory.
    """
    # List all files and directories in the specified directory
    all_items = os.listdir(directory)
    
    # Filter out only the paths of files (not directories)
    file_paths = [os.path.join(directory, item) for item in all_items if os.path.isfile(os.path.join(directory, item))]
    
    return file_paths

def get_file_paths_dict(directory):
    """
    Get a dictionary where keys are filenames and values are file paths within the directory.
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