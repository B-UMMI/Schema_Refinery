#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Purpose
-------

Accepts a matrix with results from the AlleleCall process of
chewBBACA and masks certain matrix elements by substituting
those elements with other characters.

The default masking option will substitute all ASM, ALM, NIPH,
NIPHEM, PLOT3, PLOT5, LOTSC and LNF cases with '0' and the 'INF-'
prefix of inferred alleles will always be removed to homogenize
valid allele identifiers. Passing a single word will change the
default substitution value for that word. To change specific
matrix elements, the string should be formatted as:

                    'ASM=short ALM=large LNF=pop'

which will change the default substitution value for ASM and ALM
cases, maintaining the default value of '0' to substitute remaining
cases. The 'pop' expression serves to signal that a type of case
should not be substituted.

The '=' character is used to change the default substitution value
for the right side matrix element (ASM=short --> ASM cases will be
substituted with the 'short' word) and should only be used for that
purpose.

Execution example:

    >>> python mask_matrix.py -i <allele_call_matrix_file> 
        -o <masked_output_matrix> -mc ASM=short ALM=large LNF=pop

    Will substitute ASM elements with 'short' and ALM elements with 'large'.
    LNF cases will not be substituted due to the use of the 'pop' expression.
    NIPH, NIPHEM, PLOT3, PLOT5 and LOTSC cases will be substituted by the default
    '0' character. All 'INF-' prefixes will be removed from inferred alleles.

    >>> python mask_matrix.py -i <allele_call_matrix_file>
        -o <masked_output_matrix> -mc sub

    Will change the default subtitution value from '0' to 'sub' and all 
    ASM, ALM, NIPH, NIPHEM, PLOT3, PLOT5, LOTSC and LNF cases will be 
    substituted by 'sub'.
"""


import os
import csv
import argparse


def join_iterable(iterable, delimiter='\t'):
    """
    """

    joined = delimiter.join(iterable)

    return joined


def write_text(text, output_file, mode='w'):
    """
    """

    # write matrix to output file
    with open(output_file, mode) as outfile:
        outfile.write(text+'\n')


def write_lines(lines, output_file, mode='w'):
    """ Writes a matrix to a file.

    Parameters
    ----------
    matrix_rows : list
        List of sublists where each sublist corresponds
        to one row of a AlleleCall matrix.
    output_matrix : str
        Path to the file that should be created to store
        the output matrix.

    Returns
    -------
    Writes matrix rows (each sublist of the input list is
    a row) to the output file.
    """

    # join matrix lines into chunk of text
    concat_lines = [join_iterable(line, '\t')
                    for line in lines]
    lines_text = join_iterable(concat_lines, '\n')

    write_text(lines_text, output_file, mode)


def change_masking_chars(masking_dict, masking_characters):
    """
    """

    for char in masking_characters:
        c = char.split('=')
        if c[1] != 'pop':
            masking_dict[c[0]] = c[1]
        else:
            masking_dict.pop(c[0])
        # determine if user passed a value that can be converted
        # to integer and conflict with allele ids
        try:
            new_char = int(c[1])
            if new_char != '0':
                print('WARNING: You have chosen {0} to substitute for '
                      '{1}. This might conflict with valid allele '
                      'identifiers.'.format(c[1], c[0]))
        except:
            pass

    return masking_dict


def mask_matrix(input_matrix, masking_chars_dict, output_matrix):
    """ Masks matrix values. Masked values are substituted by
        either a default character or user-defined characters.

    Parameters
    ----------
    matrix_rows : list
        List of sublists where each row has the elements of
        one row of a AlleleCall matrix.
    masking_chars_dict : dict
        Dictionary with matrix elements that should be
        substituted as keys and the characters that will
        substitute those elements as values.

    Returns
    -------
    masked_matrix : list
        List of sublists where each sublist has the elements
        of one row of the matrix, masked according to the
        character substitutions specified in the
        `masking_chars_dict` dictionary.
    """

    limit = 200
    total_masked = 0
    with open(input_matrix, 'r') as infile:
        reader = csv.reader(infile, delimiter='\t')

        masked_matrix_lines = []
        # first row is the header with loci identifiers
        masked_matrix_lines.append(next(reader))

        # for each matrix row that has an allelic profile
        exausted = False
        while exausted is False:
            try:
                row = reader.__next__()

                # remove the 'INF-' prefix from each inferred allele
                # to homogenize valid allele identifiers
                masked_row = [e.split('-')[1]
                              if 'INF-' in e
                              else e
                              for e in row]

                # substitute the matrix elements that are keys
                # in the masking dictionary
                masked_row = [masking_chars_dict[e]
                              if e in masking_chars_dict
                              else e
                              for e in masked_row]

                # append each masked row
                masked_matrix_lines.append(masked_row)
            except:
                exausted = True

            if len(masked_matrix_lines) >= limit or exausted is True:
                if len(masked_matrix_lines) > 0:
                    write_lines(masked_matrix_lines, output_matrix, 'a')
                    total_masked += len(masked_matrix_lines)
                    masked_matrix_lines = []

    return total_masked


def main(input_matrix, output_file, masking_characters):

    # determine input basename
    input_basename = os.path.basename(input_matrix)
    # remove extension that is after last '.'
    input_basename = '.'.join(input_basename.split('.')[0:-1])

    classes = ['ALM', 'ASM', 'LNF', 'NIPH',
               'NIPHEM', 'PLOT3', 'PLOT5', 'LOTSC']

    # sort list of masking characters to get character
    # to substitute default first
    if masking_characters is not None:
        masking_characters = sorted(masking_characters,
                                    key=lambda x: '=' in x)

    # define default masking character
    if masking_characters is None or '=' in masking_characters[0]:
        default_char = '0'
    elif masking_characters is not None and '=' not in masking_characters[0]:
        default_char = masking_characters[0]
        masking_characters.remove(default_char)

    masking_dict = {c: default_char for c in classes}

    # alter dictionary when users want to alter dictionary values
    # for a value different than default '0' or when cases to substitute
    # should not be all substituted by the same character
    if masking_characters is not None:
        masking_dict = change_masking_chars(masking_dict,
                                            masking_characters)
    masking_table = ['{0}: {1}'.format(k, v)
                     for k, v in masking_dict.items()]
    print('Mask table:')
    print('\n'.join(masking_table)+'\n')

    print('Masking matrix...', end='')
    # mask matrix
    masked_matrix = mask_matrix(input_matrix, masking_dict, output_file)
    print('masked matrix available at {0}'.format(os.path.abspath(output_file)))


def parse_arguments():

    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument('-i', '--input-matrix', type=str,
                        required=True, dest='input_matrix',
                        help='Path to the file with the AlleleCall '
                             'matrix (default name given by chewBBACA'
                             ' is results_alleles.tsv).')

    parser.add_argument('-o', '--output-file', type=str,
                        required=True, dest='output_file',
                        help='Path to the output file that will be '
                             'created to store the new matrix.')

    parser.add_argument('-mc', '--masking-characters', type=str,
                        nargs='+', required=False,
                        dest='masking_characters', default=None,
                        help='Define the character that will substitute'
                             ' the cases to be masked. By default, all '
                             'cases that are not valid alleles are '
                             'substituted by "0" and inferred alleles '
                             'are stripped from the "INF-" prefix. '
                             'Defining another single character will '
                             'substitute all cases with that character. '
                             'Different cases can be substituted by '
                             'different characters but the characters '
                             'must be given in the format "ASM=NewChar". '
                             'Cases that are not specified will be '
                             'substituted by the default and the value '
                             '"pop" will remove the specified cases '
                             'from the list of cases that should be '
                             'substituted.')

    args = parser.parse_args()

    return args


if __name__ == "__main__":

    args = parse_arguments()
    main(**vars(args))