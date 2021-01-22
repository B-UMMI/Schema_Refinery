#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

"""


import os
import argparse


def main(input_file, schema_dir):

    with open(input_file, 'r') as infile:
        loci_list = infile.read().splitlines()

    # list loci in schema
    schema_loci = [f for f in os.listdir(schema_dir) if '.fasta' in f]

    total = 0
    for locus in schema_loci:
        locus_id = locus.split('.fasta')[0]
        if locus_id in loci_list:
            os.remove(os.path.join(schema_dir, locus))
            os.remove(os.path.join(schema_dir, 'short', locus_id+'_short'+'.fasta'))
            total += 1

    print('Removed {0} from schema.'.format(total))


def parse_arguments():

    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument('-i', type=str, required=True,
                        dest='input_file',
                        help='')

    parser.add_argument('-s', type=str, required=True,
                        dest='schema_dir',
                        help='')

    args = parser.parse_args()

    return args


if __name__ == '__main__':

    args = parse_arguments()

    main(**vars(args))
