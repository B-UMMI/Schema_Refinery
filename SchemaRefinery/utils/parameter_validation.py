import os
import sys
import ast
import csv

try:
    from utils import constants as ct
except ModuleNotFoundError:
    from SchemaRefinery.utils import constants as ct


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


def check_minimum(value, minimum):
    """Check if a value is below a threshold value.

    Parameters
    ----------
    value : int or float
        The value to check.
    minimum : int or float
        The minimum threshold value.

    Returns
    -------
    bool
        True if the value is above or equal to the minimum, False otherwise.
    """
    return value >= minimum


def check_maximum(value, maximum):
    """Check if a value is above a threshold value.

    Parameters
    ----------
    value : int or float
        The value to check.
    maximum : int or float
        The maximum threshold value.

    Returns
    -------
    bool
        True if the value is below or equal to the maximum, False otherwise.
    """
    return value <= maximum


def check_value_interval(value, minimum, maximum):
    """Check if parameter value is contained in interval.

    Parameters
    ----------
    value : int or float
        The value to check.
    minimum : int or float
        The minimum threshold value.
    maximum : int or float
        The maximum threshold value.

    Returns
    -------
    bool
        True if the value is within the interval [minimum, maximum], False otherwise.
    """
    return check_minimum(value, minimum) and check_maximum(value, maximum)


def check_value_type(value, expected_type):
    """Check if parameter is of expected type.

    Parameters
    ----------
    value : any
        The value to check.
    expected_type : type
        The expected type of the value.

    Returns
    -------
    any
        The converted value if it matches the expected type, None otherwise.
    """
    try:
        if expected_type is bool:
            converted = tryeval(value)
        else:
            converted = expected_type(value)
        if type(converted) is expected_type:
            return converted
        else:
            return None
    except ValueError:
        return None


def check_path(value):
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


def check_in_list(values, expected_values):
    """Check if all values are in the expected values list.

    Parameters
    ----------
    values : list
        The list of values to check.
    expected_values : list
        The list of expected values.

    Returns
    -------
    list or bool
        The list of values if all are in the expected values, False otherwise.
    """
    intersection = set.intersection(set(expected_values), set(values))
    if len(intersection) < len(values):
        return False
    else:
        return values


def check_parameter(value, validate_type=None, validate_minimum=None, validate_maximum=None,
                    validate_path=False, validate_list=None):
    """Validate a value passed to a parameter.

    Parameters
    ----------
    value : any
        The value to validate.
    validate_type : type, optional
        The expected type of the value.
    validate_minimum : int or float, optional
        The minimum threshold value.
    validate_maximum : int or float, optional
        The maximum threshold value.
    validate_path : bool, optional
        Whether to check if the value is a valid path.
    validate_list : list, optional
        The list of expected values.

    Returns
    -------
    any
        The validated value if all checks pass, None otherwise.
    """
    valid = True
    if validate_type and valid:
        value = check_value_type(value, validate_type)
        if value is None:
            valid = False
    if validate_minimum and valid:
        valid_min = check_minimum(value, validate_minimum)
        if not valid_min:
            valid = False
    if validate_maximum and valid:
        valid_max = check_maximum(value, validate_maximum)
        if not valid_max:
            valid = False
    if validate_path and valid:
        valid_path = check_path(value)
        if not valid_path:
            valid = False
    if validate_list and valid:
        value = check_in_list(value.split(','), validate_list)
        if value is False:
            valid = False

    if not valid:
        return None
    else:
        return value


def validate_criteria_file(file_path, expected_criteria=ct.FILTERING_CRITERIA):
    """
    Validates modify_schema criteria input file.

    Parameter
    ---------
    file_path : str
        File path containg the criteria file.
    expected_criteria : list
        Expected Filtering criteria.

    Returns
    -------
    parameter_values : dict
        Returns dictionary containing criteria values.
    """
    
    with open(file_path, 'r', encoding='utf-8') as filters:
        criteria = dict(csv.reader(filters, delimiter='\t'))

    unexpected_keys = [x
                       for x in criteria
                       if x not in expected_criteria]

    if len(unexpected_keys) > 0:
        print("\nError: Following unexpected criteria:")
        unexpected_keys = '\n'.join(unexpected_keys)
        print(unexpected_keys)
        sys.exit()

    missing_keys = [x
                    for x in expected_criteria
                    if x not in criteria]

    if len(missing_keys) > 0:
        print("\nError: Missing following criteria:")
        missing_keys = '\n'.join(missing_keys)
        print(missing_keys)
        sys.exit()

    warnings = []
    parameter_values = {}
    for k, v in criteria.items():
        if v not in ['None', '']:
            valid = check_parameter(v, *ct.CRITERIA_ERRORS[k][1])
            if valid is not None:
                parameter_values[k] = valid
            else:
                warnings.append(ct.CRITERIA_ERRORS[k][0])
        else:
            parameter_values[k] = None

    if len(warnings) > 0:
        warning_text = '\n'.join(warnings)
        sys.exit(warning_text)
    else:
        return parameter_values
