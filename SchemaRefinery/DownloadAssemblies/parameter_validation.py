#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
"""


import os
import sys
import ast
import csv
from typing import Any, Dict, List, Optional, Union

try:
    from DownloadAssemblies import constants as ct
except ModuleNotFoundError:
    from SchemaRefinery.DownloadAssemblies import constants as ct


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