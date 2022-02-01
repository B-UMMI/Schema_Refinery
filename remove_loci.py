#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Purpose
-------
This script removes loci from a schema.

Code documentation
------------------
"""


import os
import argparse


def main(input_file, schema_directory):

    with open(input_file, 'r') as infile:
        loci_list = infile.read().splitlines()

    # list loci in schema
    schema_loci = [f
                   for f in os.listdir(schema_directory)
                   if '.fasta' in f]

    total = 0
    for locus in schema_loci:
        locus_id = locus.split('.fasta')[0]
        if locus_id in loci_list:
            # remove locus file from main schema directory
            os.remove(os.path.join(schema_directory, locus))
            # remove from short diretory
            os.remove(os.path.join(schema_directory, 'short', locus_id+'_short'+'.fasta'))
            total += 1

    print('Removed {0} loci from schema.'.format(total))


def parse_arguments():

    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument('-i', '--input-file', type=str, required=True,
                        dest='input_file',
                        help='Path to input file with the list of loci '
                             'to remove from the schema (one per line).')

    parser.add_argument('-s', '--schema-directory', type=str, required=True,
                        dest='schema_directory',
                        help='Path to the schema\'s directory.')

    args = parser.parse_args()

    return args


if __name__ == '__main__':

    args = parse_arguments()
    main(**vars(args))
