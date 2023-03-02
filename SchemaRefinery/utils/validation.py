import sys
import platform

try:
    from utils import constants as ct
except:
    from SchemaRefinery.utils import constants as ct

def validate_python_version(minimum_version=ct.MIN_PYTHON):
    """Validate Python version used to run Schema Refinery.
    Parameters
    ----------
    minimum_version : tuple
        A tuple with the Puthon version as (MAJOR, MINOR, PATCH).
        According to the rules of Semanting Versioning
        (https://semver.org/).
    Returns
    -------
    python_version : str
        Python version in format "MAJOR.MINOR.PATCH".
    Raises
    ------
    SystemExit
        - If the Python version does not meet minimum requirements
        or it was not possible to determine/detect a version.
    """
    python_version = platform.python_version()

    try:
        assert tuple(map(int, python_version.split('.'))) >= minimum_version[0]
    except AssertionError:
        print('Python version found: {} '.format(python_version))
        print('Please use version Python >= {0}'.format(minimum_version[1]))
        sys.exit(0)

    return python_version