import sys
import platform
from typing import Tuple

try:
    from utils import constants as ct
except:
    from SchemaRefinery.utils import constants as ct

def validate_python_version(minimum_version: Tuple[int, int, int] = ct.MIN_PYTHON) -> str:
    """
    Validate Python version used to run Schema Refinery.

    Parameters
    ----------
    minimum_version : Tuple[int, int, int]
        A tuple with the Python version as (MAJOR, MINOR, PATCH).
        According to the rules of Semantic Versioning
        (https://semver.org/).

    Returns
    -------
    str
        Python version in format "MAJOR.MINOR.PATCH".

    Raises
    ------
    SystemExit
        - If the Python version does not meet minimum requirements
        or it was not possible to determine/detect a version.
    """
    python_version: str = platform.python_version()

    try:
        assert tuple(map(int, python_version.split('.'))) >= minimum_version
    except AssertionError:
        print(f'Python version found: {python_version}')
        print(f'Please use Python version >= {minimum_version[0]}.{minimum_version[1]}.{minimum_version[2]}')
        sys.exit(0)

    return python_version
