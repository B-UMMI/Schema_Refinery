#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
"""


import os
import sys
import ast
import csv

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


def check_minimum(value, minimum):
    """Check if a value is below a threshold value."""
    if value < minimum:
        return False
    else:
        return True


def check_maximum(value, maximum):
    """Check if a value is above a threshold value."""
    if value > maximum:
        return False
    else:
        return True


def check_value_interval(value, minimum, maximum):
    """Check if parameter value is contained in interval."""
    if check_minimum(value) and check_maximum(value):
        return True
    else:
        False


def check_value_type(value, expected_type):
    """Check if parameter is of expected type."""
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
    """Check if a path exists."""
    if os.path.exists(value):
        return True
    else:
        return False


def check_in_list(values, expected_values):
    """"""
    intersection = set.intersection(set(expected_values), set(values))
    if len(intersection) < len(values):
        return False
    else:
        return values


def check_parameter(value, validate_type, validate_minimum, validate_maximum,
                    validate_path, validate_list):
    """Validate a value passed to a parameter."""
    valid = True
    if validate_type and valid:
        value = check_value_type(value, validate_type)
        if value is None:
            valid = False
    if validate_minimum and valid:
        valid_min = check_minimum(value, validate_minimum)
        if valid_min is None:
            valid = False
    if validate_maximum and valid:
        valid_max = check_maximum(value, validate_maximum)
        if valid_max is None:
            valid = False
    if validate_path and valid:
        valid_path = check_path(value)
        if valid_path is None:
            valid = False
    if validate_list and valid:
        value = check_in_list(value.split(','), validate_list)
        if value is None:
            valid = False

    if valid is None:
        return valid
    else:
        return value


def validate_criteria_file(file_path, expected_criteria=ct.FILTERING_CRITERIA):
    """"""
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

def validate_modify_schema_input(file_path, expected_commands=ct.ACCEPTED_COMMANDS):
    input_list = []
    with open(file_path, 'r', newline='') as file:
        reader = csv.reader(file, delimiter='\t')
        for row in reader:
            non_empty_row = [value for value in row if value]
            if non_empty_row:
                input_list.append(non_empty_row)

    unexpected_keys = [x
                       for x in input_list
                       if x[0] not in expected_commands]

    if len(unexpected_keys) > 0:
        print("\nError: Following unexpected commands:")
        unexpected_keys = '\n'.join(unexpected_keys)
        print(unexpected_keys)
        sys.exit()

    warnings = []
    parameter_values = []
    for input_line in input_list:
        command = input_line[0].lower()
        parameters = input_line[1:]
        validated_values = []
        if command in ['merge', 'remove']:
            for parameter in parameters:
                validated_values.append(check_parameter(parameter, *ct.INPUT_ERRORS[command][1]))
        else:
            if (len(parameters) % 2) != 0:
                for new_loci_name, interval in zip(*[iter(parameter)]*2):
                    validated_values.append(check_parameter(new_loci_name, *ct.INPUT_ERRORS[command + "_id"][1]))

                    try:
                        valid = []
                        for i in interval.split('-'):
                            valid.append(check_parameter(i, *ct.INPUT_ERRORS[command + "_size"][1]))
                        
                        if all(valid is not None):
                            validated_values.append(interval)
                    except:
                        if ct.SPLIT_MISSING_MINUS not in warnings:
                            warnings.append(ct.SPLIT_MISSING_MINUS)

            elif ct.SPLIT_MUST_BE_EVEN not in warnings:
                warnings.append(ct.SPLIT_MUST_BE_EVEN)
        if all(validated_values):
            parameter_values.append([command] + validated_values)
        elif ct.CRITERIA_ERRORS[command][0] not in warnings:
            warnings.append(ct.CRITERIA_ERRORS[command][0])

    if len(warnings) > 0:
        warning_text = '\n'.join(warnings)
        sys.exit(warning_text)
    else:
        return parameter_values


        

