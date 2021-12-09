#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Purpose
-------
This script determines if any locus in a schema does
not have sequential allele identifiers.

Code documentation
------------------
"""


import os
import argparse

from Bio import SeqIO


def main(schema_directory):

    schema_loci = [os.path.join(schema_directory, file)
                   for file in os.listdir(schema_directory)
                   if '.fasta' in file]

    missing_ids = {}
    for locus in schema_loci:
        locus_id = os.path.basename(locus).split('.fasta')[0]
        allele_ids = [int((rec.id).split('_')[-1])
                      for rec in SeqIO.parse(locus, 'fasta')]

        # determine if there are missing ids
        min_id = 1
        max_id = max(allele_ids)
        valid_range = list(range(min_id, max_id+1))
        missing = sorted(set(valid_range).difference(allele_ids))

        if len(missing) > 0:
            missing_ids[locus_id] = missing
            print(locus, missing)


def parse_arguments():

    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument('-s', type=str, required=True,
                        dest='schema_directory',
                        help='Path to the schema\'s directory.')

    args = parser.parse_args()

    return args


if __name__ == '__main__':

    args = parse_arguments()
    main(**vars(args))
