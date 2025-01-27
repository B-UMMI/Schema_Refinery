import shutil
from typing import Optional

from SchemaRefinery.utils import print_functions as pf

def get_tool_path(name: str) -> str:
    """
    Get the path of the specified tool.

    Parameters
    ----------
    name : str
        The name of the tool.

    Returns
    -------
    Optional[str]
        The path of the tool if it is found in the PATH. If the tool is not found, returns None.
    
    Raises
    ------
    FileNotFoundError
        If the tool is not found in the PATH.
    """
    
    path: Optional[str] = shutil.which(name)

    if path is None:
        pf.print_message(f"Error: {name} not found in PATH.", "error")
        raise FileNotFoundError
    else:
        return path