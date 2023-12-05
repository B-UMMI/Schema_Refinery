import os

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