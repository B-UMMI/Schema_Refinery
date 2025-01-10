import os
import sys
import ast
import csv
import argparse
import platform
from typing import Tuple, Dict, Any, List, Union, Optional

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

def tryeval(val):
    """
    Evaluates the type of the input.

    Parameter
    ---------
    val : any type

    Returns
    -------
    val : any type
        converted to the right type.
    """

    try:
        val = ast.literal_eval(val)
    except ValueError:
        pass
    return val


def check_minimum(value: Union[int, float], minimum: Union[int, float]) -> bool:
    """Check if a value is below a threshold value.

    Parameters
    ----------
    value : Union[int, float]
        The value to check.
    minimum : Union[int, float]
        The minimum threshold value.

    Returns
    -------
    bool
        True if value is above or equal to minimum, False otherwise.
    """
    return value >= minimum


def check_maximum(value: Union[int, float], maximum: Union[int, float]) -> bool:
    """Check if a value is above a threshold value.

    Parameters
    ----------
    value : Union[int, float]
        The value to check.
    maximum : Union[int, float]
        The maximum threshold value.

    Returns
    -------
    bool
        True if value is below or equal to maximum, False otherwise.
    """
    return value <= maximum


def check_value_interval(value: Union[int, float], minimum: Union[int, float], maximum: Union[int, float]) -> bool:
    """Check if parameter value is contained in interval.

    Parameters
    ----------
    value : Union[int, float]
        The value to check.
    minimum : Union[int, float]
        The minimum threshold value.
    maximum : Union[int, float]
        The maximum threshold value.

    Returns
    -------
    bool
        True if value is within the interval, False otherwise.
    """
    return check_minimum(value, minimum) and check_maximum(value, maximum)


def check_value_type(value: Any, expected_type: type) -> Optional[Any]:
    """Check if parameter is of expected type.

    Parameters
    ----------
    value : Any
        The value to check.
    expected_type : type
        The expected type of the value.

    Returns
    -------
    Optional[Any]
        The converted value if it matches the expected type, None otherwise.
    """
    try:
        if expected_type is bool:
            converted = tryeval(value)
        else:
            converted = expected_type(value)
        if isinstance(converted, expected_type):
            return converted
        else:
            return None
    except ValueError:
        return None


def check_path(value: str) -> bool:
    """Check if a path exists.

    Parameters
    ----------
    value : str
        The path to check.

    Returns
    -------
    bool
        True if the path exists, False otherwise.
    """
    return os.path.exists(value)


def check_in_list(values: List[str], expected_values: List[str]) -> Union[bool, List[str]]:
    """Check if all values are in the expected values list.

    Parameters
    ----------
    values : List[str]
        The values to check.
    expected_values : List[str]
        The list of expected values.

    Returns
    -------
    Union[bool, List[str]]
        The values if all are in the expected values list, False otherwise.
    """
    intersection = set.intersection(set(expected_values), set(values))
    if len(intersection) < len(values):
        return False
    else:
        return values


def check_parameter(value: Any, validate_type: Optional[type], validate_minimum: Optional[Union[int, float]],
                    validate_maximum: Optional[Union[int, float]], validate_path: bool, validate_list: Optional[List[str]]) -> Union[bool, Any]:
    """Validate a value passed to a parameter.

    Parameters
    ----------
    value : Any
        The value to validate.
    validate_type : Optional[type]
        The expected type of the value.
    validate_minimum : Optional[Union[int, float]]
        The minimum threshold value.
    validate_maximum : Optional[Union[int, float]]
        The maximum threshold value.
    validate_path : bool
        Whether to validate the value as a path.
    validate_list : Optional[List[str]]
        The list of expected values.

    Returns
    -------
    Union[bool, Any]
        The validated value if valid, False otherwise.
    """
    valid: bool = True
    if validate_type and valid:
        value = check_value_type(value, validate_type)
        if value is None:
            valid = False
    if validate_minimum and valid:
        valid_min: bool = check_minimum(value, validate_minimum)
        if not valid_min:
            valid = False
    if validate_maximum and valid:
        valid_max: bool = check_maximum(value, validate_maximum)
        if not valid_max:
            valid = False
    if validate_path and valid:
        valid_path: bool = check_path(value)
        if not valid_path:
            valid = False
    if validate_list and valid:
        value = check_in_list(value.split(','), validate_list)
        if value is None:
            valid = False

    return value if valid else False


def validate_criteria_file(file_path: str, expected_criteria: Dict[str, Any] = ct.FILTERING_CRITERIA) -> Dict[str, Any]:
    """Validates initial input criteria arguments file to be according to the desired format.

    Parameters
    ----------
    file_path : str
        File path to the criteria file.
    expected_criteria : Dict[str, Any]
        Contains the type and format that criteria file is supposed to have.

    Returns
    -------
    Dict[str, Any]
        Dictionary that contains the parameters values extracted from criteria file.
    """
    with open(file_path, 'r', encoding='utf-8') as filters:
        criteria: Dict[str, str] = dict(csv.reader(filters, delimiter='\t'))

    unexpected_keys: List[str] = [x for x in criteria if x not in expected_criteria]

    if unexpected_keys:
        print("\nError: Following unexpected criteria:")
        print('\n'.join(unexpected_keys))
        sys.exit()

    missing_keys: List[str] = [x for x in expected_criteria if x not in criteria]

    if missing_keys:
        print("\nError: Missing following criteria:")
        print('\n'.join(missing_keys))
        sys.exit()

    warnings: List[str] = []
    parameter_values: Dict[str, Any] = {}
    for k, v in criteria.items():
        if v not in ['None', '']:
            valid: Union[bool, Any] = check_parameter(v, *ct.CRITERIA_ERRORS[k][1])
            if valid is not None:
                parameter_values[k] = valid
            else:
                warnings.append(ct.CRITERIA_ERRORS[k][0])
        else:
            parameter_values[k] = None

    if warnings:
        sys.exit('\n'.join(warnings))
    else:
        return parameter_values


def validate_schema_annotation_module_arguments(args: argparse.Namespace) -> None:
    """
    Validate the arguments passed to the schema annotation module.

    Raises
    ------
    SystemExit
        - If the arguments are invalid.
    """
    # Verify if files or directories exist
    verify_schema_sctructure(args.schema_directory)

    # Arguments to match uniprot-proteomes
    if 'uniprot-proteomes' in args.annotation_options:
        if args.proteome_ids and not args.proteome_table:
            sys.exit("\nError: 'proteome_ids' can only be used with '--proteome-table' and 'uniprot-proteomes' annotation option.")
        if not args.proteome_table:
            sys.exit("\nError: 'proteome-table' is required with 'uniprot-proteomes' annotation option.")
        # Verify if files or directories exist
        if args.proteome_table:
            # Verify if files or directories exist
            verify_path_exists(args.proteome_table, 'file')
    else:
        if any([args.proteome_ids, args.proteome_table]):
            sys.exit("\nError: 'proteome_ids' and 'proteome-table' can only be used with '--annotation-options uniprot-proteomes'.")

    # Arguments to match genbank
    if 'genbank' in args.annotation_options:
        if args.genbank_ids_to_add and not args.genbank_files:
            sys.exit("\nError: 'genbank-ids-to-add' can only be used with '--genbank-files' and 'genbank' annotation option.")
        if not args.genbank_files:
            sys.exit("\nError: 'genbank-files' is required with 'genbank' annotation option.")
        # Verify if files or directories exist
        if args.genbank_files:
            # Verify if files or directories exist
            verify_path_exists(args.genbank_files, 'directory')
            # Verify if the GenBank files directory is empty
            if_genbank_files_empty = not os.listdir(args.genbank_files)
            if if_genbank_files_empty:
                sys.exit("\nError: The GenBank files directory is empty.")
    else:
        if any([args.genbank_files, args.genbank_ids_to_add]):
            sys.exit("\nError: 'genbank-files' and 'genbank-ids-to-add' can only be used with '--annotation-options genbank'.")

    # Arguments to match schemas
    if 'match-schemas' in args.annotation_options:
        if any([args.subject_schema, args.subject_annotations, args.processing_mode]) and not all([args.subject_schema, args.subject_annotations, args.processing_mode]):
            missing_args = []
            if not args.subject_schema:
                missing_args.append('subject_schema')
            if not args.subject_annotations:
                missing_args.append('subject_annotations')
            if not args.processing_mode:
                missing_args.append('processing_mode')

            sys.exit(f"\nError: Missing required arguments: {', '.join(missing_args)}. All of 'subject_schema', 'subject_annotations', and 'processing_mode' must be provided together.")
        
        if args.subject_schema:
            verify_path_exists(args.subject_schema, 'directory')

        if args.subject_annotations:
            verify_path_exists(args.subject_annotations, 'file')
    else:
        if any([args.subject_schema, args.subject_annotations, args.processing_mode]):
            sys.exit("\nError: 'subject_schema', 'subject_annotations', and 'processing_mode' can only be used with '--annotation-options match-schemas'.")


def validate_download_assemblies_module_arguments(args: argparse.Namespace) -> None:
    """
    Validate the arguments passed to the download assemblies module.

    Raises
    ------
    SystemExit
        - If the arguments are invalid
    """
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
