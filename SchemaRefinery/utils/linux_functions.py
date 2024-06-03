import shutil

def get_tool_path(name):
    """
    Get the path of the specified tool.

    Parameters
    ----------
    name : str
        The name of the tool.

    Returns
    -------
    return : str or None
        The path of the tool if it is found in the PATH. If the tool is not found, returns None.
    """

    # Get the path of the tool
    return shutil.which(name)