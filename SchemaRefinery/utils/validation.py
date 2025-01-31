import os
import sys
import ast
import csv
import argparse
import platform
from typing import Tuple, Dict, Any, List, Union, Optional

try:
    from utils import (constants as ct,
                       print_functions as pf)
except:
    from SchemaRefinery.utils import (constants as ct,
                                      print_functions as pf)

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
        pf.print_message(f'Python version found: {python_version}', 'info')
        pf.print_message(f'Please use Python version >= {minimum_version[0]}.{minimum_version[1]}.{minimum_version[2]}', 'error')
        sys.exit(0)

    return python_version

def validate_system_max_cpus_number(given_number_of_cpus: int) -> int:
    """
    Validate the number of CPUs available in the system and return the optimal number of CPUs to use.

    Parameters
    ----------
    given_number_of_cpus : int
        The number of CPUs to use.

    Returns
    -------
    int
        The number of CPUs to use.
    """
    max_cpus: Optional[int] = os.cpu_count()
    if max_cpus and max_cpus < given_number_of_cpus:
        optimal_cpus: int = max_cpus - 1
        print(f'The number of CPUs available in the system is less than the given number of CPUs: {max_cpus} < {given_number_of_cpus}, setting cpu count to {optimal_cpus}', 'warning')
        return optimal_cpus
    return given_number_of_cpus

def verify_path_exists(path: str, path_type: str, errors: List[str]) -> None:
    """
    Verify if a file or directory exists and append errors to the provided list.

    Parameters
    ----------
    path : str
        The path to the file or directory.
    path_type : str
        The type of path ('file' or 'directory').
    errors : List[str]
        The list to append error messages to.
    """
    if not os.path.exists(path):
        errors.append(f"The specified {path_type} does not exist: {path}")

def verify_schema_structure(schema_directory: str, errors: List[str]) -> None:
    """
    Verify the structure of the schema directory and append errors to the provided list.

    Parameters
    ----------
    schema_directory : str
        The path to the schema directory.
    errors : List[str]
        The list to append error messages to.
    """
    if not os.path.exists(schema_directory):
        errors.append(f"The schema directory does not exist: {schema_directory}")

    files = [os.path.join(schema_directory, f) for f in os.listdir(schema_directory) if os.path.isfile(os.path.join(schema_directory, f))]
    folders = [os.path.join(schema_directory, f) for f in os.listdir(schema_directory) if os.path.isdir(os.path.join(schema_directory, f))]
    
    if len(files) == 0:
        errors.append(f"The schema directory is empty: {schema_directory}")
    
    if len(folders) == 0:
        errors.append(f"The schema directory is missing the schema short: {schema_directory}")
    
    short_schema_directory = os.path.join(schema_directory, 'short')
    if not os.path.exists(short_schema_directory):
        errors.append(f"The schema short directory does not exist: {short_schema_directory}")

    short_files = [os.path.join(short_schema_directory, f) for f in os.listdir(short_schema_directory) if os.path.isfile(os.path.join(short_schema_directory, f))]

    if len(short_files) == 0:
        errors.append(f"The schema short directory is empty: {short_schema_directory}")

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
        pf.print_message("Following unexpected criteria:", "error")
        unexpected_keys_string: str = '\n'.join(unexpected_keys)
        pf.print_message(unexpected_keys_string, None)
        sys.exit()

    missing_keys: List[str] = [x for x in expected_criteria if x not in criteria]

    if missing_keys:
        pf.print_message("Missing following criteria:", "error")
        missing_keys_string: str = '\n'.join(missing_keys)
        pf.print_message(missing_keys_string, None)
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
        pf.print_message("The following errors were found:", "error")
        pf.print_message('\n'.join(warnings), None)
        sys.exit()
    else:
        return parameter_values

def validate_download_assemblies_module_arguments(args: argparse.Namespace) -> None:
    """
    Validate the arguments passed to the download assemblies module.

    Parameters
    ----------
    args : argparse.Namespace
        The arguments passed to the download assemblies module.

    Raises
    ------
    SystemExit
        - If the arguments are invalid
    """
    errors: List[str] = []

    # Verify if files or directories exist
    verify_path_exists(args.output_directory, 'directory', errors)

    # Ensure that ENA661K is not used with an input table
    if args.input_table is not None and 'ENA661K' in args.database:
        errors.append("\nError: Only assemblies from NCBI can be fetched from an input file. ENA661K was parsed.")
    
    if args.filtering_criteria.taxon and args.filtering_criteria.input_table:
        errors.append("\nError: Taxon and input table are mutually exclusive.")

    if not args.filtering_criteria.taxon and not args.filtering_criteria.input_table:
        errors.append("\nError: Must provide either taxon or input table.")
    
    if args.filtering_criteria.input_table and 'NCBI' not in args.database: 
        errors.append("\nError: Input table can only be used with NCBI database.")

    elif args.filtering_criteria.input_table:
        verify_path_exists(args.filtering_criteria.input_table, 'file', errors)
    
    if args.threads <= 0:
        errors.append("\nError: 'threads' must be a value greater than 0.")
    
    if args.retry <= 0:
        errors.append("\nError: 'retry' must be a value greater than 0.")

    # Display all errors at once if there are any
    if errors:
        pf.print_message("The following errors were found:", "error")
        pf.print_message("\n".join(errors))
        sys.exit()

def validate_schema_annotation_module_arguments(args: argparse.Namespace) -> None:
    """
    Validate the arguments passed to the schema annotation module.

    Parameters
    ----------
    args : argparse.Namespace
        The arguments passed to the schema annotation module.

    Raises
    ------
    SystemExit
        - If the arguments are invalid.
    """
    errors: List[str] = []

    # Verify if files or directories exist
    verify_schema_structure(args.schema_directory, errors)
    verify_path_exists(args.output_directory, 'directory', errors)

    # Chewie annotations
    if args.chewie_annotations:
        verify_path_exists(args.chewie_annotations, 'file', errors)
    
    if args.bsr <= 0 or args.bsr >= 1:
        errors.append("\nError: 'bsr' must be a value between 0 and 1.")
    
    if args.threads <= 0:
        errors.append("\nError: 'threads' must be a value greater than 0.")
    
    if args.cpu <= 0:
        errors.append("\nError: 'cpu' must be a value greater than 0.")
    else:
        args.cpu = validate_system_max_cpus_number(args.cpu)

    if args.retry <= 0:
        errors.append("\nError: 'retry' must be a value greater than 0.")
    
    if args.translation_table < 0 or args.translation_table >= 25:
        errors.append("\nError: 'translation-table' must be a value between 0 and 25.")

    if args.clustering_sim <= 0 or args.clustering_sim >= 1:
        errors.append("\nError: 'clustering-sim' must be a value between 0 and 1.")
    
    if args.clustering_cov <= 0 or args.clustering_cov >= 1:
        errors.append("\nError: 'clustering-cov' must be a value between 0 and 1.")
    
    if args.size_ratio <= 0 or args.size_ratio >= 1:
        errors.append("\nError: 'size-ratio' must be a value between 0 and 1.")

    # Arguments to match uniprot-proteomes
    if 'uniprot-proteomes' in args.annotation_options:
        if args.proteome_ids_to_add and not args.proteome_table:
            errors.append("\nError: 'proteome-ids' can only be used with '--proteome-table' and 'uniprot-proteomes' annotation option.")
        if not args.proteome_table:
            errors.append("\nError: 'proteome-table' is required with 'uniprot-proteomes' annotation option.")
        # Verify if files or directories exist
        if args.proteome_table:
            verify_path_exists(args.proteome_table, 'file', errors)
    else:
        if any([args.proteome_ids_to_add, args.proteome_table]):
            errors.append("\nError: 'proteome-ids' and 'proteome-table' can only be used with '--annotation-options uniprot-proteomes'.")

    # Arguments to match genbank
    if 'genbank' in args.annotation_options:
        if args.genbank_ids_to_add and not args.genbank_files:
            errors.append("\nError: 'genbank-ids-to-add' can only be used with '--genbank-files' and 'genbank' annotation option.")
        if not args.genbank_files:
            errors.append("\nError: 'genbank-files' is required with 'genbank' annotation option.")
        # Verify if files or directories exist
        if args.genbank_files:
            verify_path_exists(args.genbank_files, 'directory', errors)
            # Verify if the GenBank files directory is empty
            if_genbank_files_empty = not os.listdir(args.genbank_files)
            if if_genbank_files_empty:
                errors.append("\nError: The GenBank files directory is empty.")
    else:
        if any([args.genbank_files, args.genbank_ids_to_add]):
            errors.append("\nError: 'genbank-files' and 'genbank-ids-to-add' can only be used with '--annotation-options genbank'.")

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

            errors.append(f"\nError: Missing required arguments: {', '.join(missing_args)}. All of 'subject-schema', 'subject-annotations', and 'processing-mode' must be provided together.")
        
        if args.subject_schema:
            verify_path_exists(args.subject_schema, 'directory', errors)

        if args.subject_annotations:
            verify_path_exists(args.subject_annotations, 'file', errors)
        # Verify if best-annotations-bsr is a value between 0 and 1
        if args.best_annotations_bsr <= 0 or args.best_annotations_bsr > 1:
            errors.append("\nError: 'best-annotations-bsr' must be a value between 0 and 1.")
    else:
        if any([args.subject_schema, args.subject_annotations, args.processing_mode]):
            errors.append("\nError: 'subject-schema', 'subject-annotations', and 'processing-mode' can only be used with '--annotation-options match-schemas'.")

    # Display all errors at once if there are any
    if errors:
        pf.print_message("The following errors were found:", "error")
        pf.print_message("\n".join(errors))
        sys.exit()

def validate_identify_spurious_genes_module_arguments(args: argparse.Namespace) -> None:
    """
    Validate the arguments passed to the identify spurious genes module.

    Parameters
    ----------
    args : argparse.Namespace
        The arguments passed to the identify spurious genes module.

    Raises
    ------
    SystemExit
        - If the arguments are invalid
    """

    errors: List[str] = []

    # Verify if files or directories exist
    verify_schema_structure(args.schema_directory, errors)
    verify_path_exists(args.output_directory, 'directory', errors)
    verify_path_exists(args.allelecall_directory, 'file', errors)

    if args.possible_new_loci:
        verify_schema_structure(args.possible_new_loci, errors)

    if args.run_mode == 'unclassified-cds':
        if args.possible_new_loci:
            errors.append("\nError: 'possible-new-loci' cannot be used with 'unclassified_cds' run mode.")
    
    if args.alignment_ratio_threshold < 0 or args.alignment_ratio_threshold >= 1:
        errors.append("\nError: 'alignment-ratio-threshold' must be a value between 0 and 1.")
    
    if args.pident_threshold < 0 or args.pident_threshold >= 100:
        errors.append("\nError: 'pident-threshold' must be a value between 0 and 100.")
    
    if args.clustering_sim_threshold < 0 or args.clustering_sim_threshold >= 1:
        errors.append("\nError: 'clustering-sim-threshold' must be a value between 0 and 1.")
    
    if args.clustering_cov_threshold < 0 or args.clustering_cov_threshold >= 1:
        errors.append("\nError: 'clustering-cov-threshold' must be a value between 0 and 1.")

    if args.genome_presence:
        if args.genome_presence < 0:
            errors.append("\nError: 'genome-presence' must be a value greater than 0.")
    
    if args.absolute_size < 0:
        errors.append("\nError: 'absolute-size' must be a value greater than 0.")
    
    if args.translation_table < 0 or args.translation_table >= 25:
        errors.append("\nError: 'translation-table' must be a value between 0 and 25.")
    
    if args.bsr < 0 or args.bsr >= 1:
        errors.append("\nError: 'bsr' must be a value between 0 and 1.")
    
    if args.size_ratio < 0 or args.size_ratio >= 1:
        errors.append("\nError: 'size-ratio' must be a value between 0 and 1.")
    
    if args.cpu <= 0:
        errors.append("\nError: 'cpu' must be a value greater than 0.")
    else:
        args.cpu = validate_system_max_cpus_number(args.cpu)

    # Display all errors at once if there are any
    if errors:
        pf.print_message("The following errors were found:", "error")
        pf.print_message("\n".join(errors))
        sys.exit()

def validate_adapt_loci_module_arguments(args: argparse.Namespace) -> None:
    """
    Validate the arguments passed to the adapt loci module.

    Parameters
    ----------
    args : argparse.Namespace
        The arguments passed to the adapt loci module.

    Raises
    ------
    SystemExit
        - If the arguments are invalid
    """

    errors: List[str] = []

    # Verify if files or directories exist
    verify_path_exists(args.input_file, 'file', errors)
    verify_path_exists(args.output_directory, 'directory', errors)

    if args.cpu <= 0:
        errors.append("Error: 'cpu' must be a value greater than 0.")
    else:
        args.cpu = validate_system_max_cpus_number(args.cpu)

    if args.bsr < 0 or args.bsr >= 1:
        errors.append("Error: 'bsr' must be a value between 0 and 1.")

    if args.translation_table < 0 or args.translation_table >= 25:
        errors.append("Error: 'translation-table' must be a value between 0 and 25.")

    # Display all errors at once if there are any
    if errors:
        sys.exit("\n".join(errors))
    
def validate_identify_paralogous_loci_arguments(args: argparse.Namespace) -> None:
    """
    Validate the arguments passed to the identify paralogous loci module.

    Parameters
    ----------
    args : argparse.Namespace
        The arguments passed to the identify paralogous loci module.

    Raises
    ------
    SystemExit
        - If the arguments are invalid
    """

    errors: List[str] = []

    # Verify if files or directories exist
    verify_schema_structure(args.schema_directory, errors)
    verify_path_exists(args.output_directory, 'directory', errors)

    if args.cpu <= 0:
        errors.append("Error: 'cpu' must be a value greater than 0.")
    else:
        args.cpu = validate_system_max_cpus_number(args.cpu)

    if args.bsr < 0 or args.bsr >= 1:
        errors.append("Error: 'bsr' must be a value between 0 and 1.")

    if args.translation_table < 0 or args.translation_table >= 25:
        errors.append("Error: 'translation-table' must be a value between 0 and 25.")
    
    if args.size_threshold < 0 or args.size_threshold >= 1:
        errors.append("Error: 'size-threshold' must be a value between 0 and 1.")

    # Display all errors at once if there are any
    if errors:
        pf.print_message("The following errors were found:", "error")
        pf.print_message("\n".join(errors))
        sys.exit()

def validate_match_schemas(args: argparse.Namespace) -> None:
    """
    Validate the arguments passed to the match schemas module.

    Parameters
    ----------
    args : argparse.Namespace
        The arguments passed to the match schemas module.

    Raises
    ------
    SystemExit
        - If the arguments are invalid
    """
    errors: List[str] = []

    # Verify if files or directories exist
    verify_schema_structure(args.query_schema_directory, errors)

    verify_schema_structure(args.subject_schema_directory, errors)

    verify_path_exists(args.output_directory, 'directory', errors)

    if args.cpu <= 0:
        errors.append("Error: 'cpu' must be a value greater than 0.")
    else:
        args.cpu = validate_system_max_cpus_number(args.cpu)

    if args.bsr < 0 or args.bsr >= 1:
        errors.append("Error: 'bsr' must be a value between 0 and 1.")

    if args.translation_table < 0 or args.translation_table >= 25:
        errors.append("Error: 'translation-table' must be a value between 0 and 25.")

    # Display all errors at once if there are any
    if errors:
        pf.print_message("The following errors were found:", "error")
        pf.print_message("\n".join(errors))
        sys.exit()

def validate_create_schema_structure(args: argparse.Namespace) -> None:
    """
    Validate the arguments passed to the create schema structure module.

    Parameters
    ----------
    args : argparse.Namespace
        The arguments passed to the create schema structure module.

    Raises
    ------
    SystemExit
        - If the arguments are invalid
    """
    errors: List[str] = []

    # Verify if files or directories exist
    verify_path_exists(args.recommendations_file, 'file', errors)
    verify_path_exists(args.fastas_folder, 'directory', errors)
    verify_path_exists(args.output_directory, 'directory', errors)

    if args.cpu <= 0:
        errors.append("Error: 'cpu' must be a value greater than 0.")
    else:
        args.cpu = validate_system_max_cpus_number(args.cpu)
    
    if args.bsr < 0 or args.bsr >= 1:
        errors.append("Error: 'bsr' must be a value between 0 and 1.")

    if args.translation_table < 0 or args.translation_table >= 25:
        errors.append("Error: 'translation-table' must be a value between 0 and 25.")
    
    # Display all errors at once if there are any
    if errors:
        pf.print_message("The following errors were found:", "error")
        pf.print_message("\n".join(errors))
        sys.exit()