#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Purpose
-------
This script determines if a locus has duplicated alleles and
if different loci have alleles in common.

Code documentation
------------------
"""


import os
import hashlib
import argparse

from Bio import SeqIO


def main(schema_directory):

    # get list of loci in schema
    schema_loci = [os.path.join(schema_directory, file)
                   for file in os.listdir(schema_directory)
                   if '.fasta' in file]

    loci_hashes = {}
    for locus in schema_loci:
        locus_id = (os.path.basename(locus)).split('.fasta')[0]
        record_generator = SeqIO.parse(locus, 'fasta')
        locus_hashes = {}
        for record in record_generator:
            allele_id = (record.id).split('_')[-1]
            sequence = str(record.seq)
            seq_hash = hashlib.sha256(sequence.encode('utf-8')).hexdigest()
            locus_hashes.setdefault(seq_hash, []).append((locus_id, allele_id))
            loci_hashes.setdefault(seq_hash, []).append((locus_id, allele_id))

        # determine if locus has duplicated alleles
        locus_duplicates = [(k, v)
                            for k, v in locus_hashes.items()
                            if len(v) > 1]
        if len(locus_duplicates) > 0:
            print('Locus {0} has duplicates:\n{1}'.format(locus_id, locus_duplicates))

    # determine if different loci have alleles in common
    loci_duplicates = []
    i = 0
    for k, v in loci_hashes.items():
        distinct_loci = set([l[0] for l in v])
        if len(distinct_loci) > 1:
            print('Loci {0} have alleles in common:\n{1}'.format(distinct_loci, v))
            loci_duplicates.append(v)
        i += 1


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
