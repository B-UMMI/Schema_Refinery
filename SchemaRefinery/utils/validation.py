import os
import sys
import platform
from typing import Tuple, Dict, Any

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


def verify_path_exists(path: str, path_type: str) -> None:
    """
    Verify if a file or directory exists.

    Parameters
    ----------
    path : str
        The path to the file or directory.
    path_type : str
        The type of path ('file' or 'directory').

    Raises
    ------
    FileNotFoundError
        If the file or directory does not exist.
    """
    if not os.path.exists(path):
        raise FileNotFoundError(f"The specified {path_type} does not exist: {path}")

def verify_schema_sctructure(schema_directory: str) -> None:
    """
    Verify if the schema directory has the correct structure.

    Parameters
    ----------
    schema_directory : str
        The path to the schema directory.

    Raises
    ------
    FileNotFoundError
        If the schema directory does not have the correct structure.
    """
    if not os.path.exists(os.path.join(schema_directory)):
        raise FileNotFoundError(f"The schema directory does not exist: {schema_directory}")

    files = [os.path.join(schema_directory, f) for f in os.listdir(schema_directory) if os.path.isfile(os.path.join(schema_directory, f))]
    folders = [os.path.join(schema_directory, f) for f in os.listdir(schema_directory) if os.path.isdir(os.path.join(schema_directory, f))]
    if len(files) == 0:
        raise FileNotFoundError(f"The schema directory is empty: {schema_directory}")
    
    if len(folders) == 0:
        raise FileNotFoundError(f"The schema directory is missing the schema short: {schema_directory}")
    
    short_schema_directory = os.path.join(schema_directory, 'short')
    short_files = [os.path.join(short_schema_directory, f) for f in os.listdir(short_schema_directory) if os.path.isfile(os.path.join(short_schema_directory, f))]

    if len(short_files) == 0:
        raise FileNotFoundError(f"The schema short directory is empty: {short_schema_directory}")

def validate_schema_annotation_module_arguments(args: dict) -> None:
    """
    Validate the arguments passed to the schema annotation module.

    Raises
    ------
    SystemExit
        - If the arguments are invalid.
    """
    # Verify if files or directories exist
    verify_schema_sctructure(args.schema_directory)

    # Verify if files or directories exist
    if args.proteome_table:
        verify_path_exists(args.proteome_table, 'file')

    # Verify if files or directories exist
    if args.genbank_files:
        # Verify if files or directories exist
        verify_path_exists(args.genbank_files, 'directory')
        # Verify if the GenBank files directory is empty
        if_genbank_files_empty = not os.listdir(args.genbank_files)
        if if_genbank_files_empty:
            sys.exit("\nError: The GenBank files directory is empty.")

    # Verify if files or directories exist
    if args.chewie_annotations:
        for annotation in args.chewie_annotations:
            verify_path_exists(annotation, 'file')

    if args.subject_schema:
        verify_path_exists(args.subject_schema, 'directory')

    # Check for mutually inclusive options
    if not args.processing_mode and args.subject_schema:
        verify_path_exists(args.subject_schema, 'directory')
        sys.exit("-pm --processing-mode is required when you want to add match with subject schema.")

    # Check for mutually inclusive options
    if args.processing_mode and not args.subject_schema:
        sys.exit("-ss --subject_schema is required when you want to add processing mode.")

    # Check for mutually inclusive options
    if not args.extra_genbank_table_columns and args.genbank_files:
        sys.exit("-ss --subject_schema is required when you want to add processing mode.")

def validate_download_assemblies_module_arguments(args: dict) -> None:
    """
    Validate the arguments passed to the download assemblies module.

    Raises
    ------
    SystemExit
        - If the arguments are invalid
    """
    verify_schema_sctructure(args.schema_directory)
    # Verify if files or directories exist
    if args.input_table:
        verify_path_exists(args.input_table, 'file')

    # Check for mutually exclusive options
    if args.input_table is not None and args.taxon is not None:
        sys.exit("\nError: Downloading from input table or downloading by taxon are mutually exclusive.")

    # Ensure that either input table or taxon name is provided
    if args.input_table is None and args.taxon is None:
        sys.exit("\nError: Must provide an input table or a taxon name.")

    # Ensure that ENA661K is not used with an input table
    if args.input_table is not None and 'ENA661K' in args.database:
        sys.exit("\nError: Only assemblies from NCBI can be fetched from an input file. ENA661K was parsed.")

    # Handle filtering criteria
    if args.filtering_criteria:
        criteria: Dict[str, Any] = args.filtering_criteria
        if criteria['assembly_source'] == ['GenBank'] and (criteria['verify_status'] is True or criteria['verify_status'] is None):
            sys.exit("\nError: Assembly status can only be verified for assemblies obtained from RefSeq (Set to False, Default(None) = True)")
    else:
        print("\nNo filtering criteria provided.")