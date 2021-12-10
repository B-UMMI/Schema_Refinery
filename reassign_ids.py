#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Purpose
-------
This script reassigns allele identifiers to ensure that
the set of alleles for a locus have sequential identifiers.

Code documentation
------------------
"""


import argparse

from Bio import SeqIO


def main(input_file, output_file, start_id):

    # import sequences
    sequence_generator = SeqIO.parse(input_file, 'fasta')

    start = start_id
    new_records = []
    for record in sequence_generator:
        recid = record.id
        sequence = str(record.seq)
        new_recid = '_'.join(recid.split('_')[0:-1]) + '_' + str(start)
        new_record = '>{0}\n{1}'.format(new_recid, sequence)
        new_records.append(new_record)
        start += 1

    records_text = '\n'.join(new_records)

    with open(output_file, 'w') as outfile:
        outfile.write(records_text+'\n')


def parse_arguments():

    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument('-i', '--input-file', type=str,
                        required=True, dest='input_file',
                        help='Path to input Fasta file.')

    parser.add_argument('-o', '--output-file', type=str,
                        required=True, dest='output_file',
                        help='Path to output file.')

    parser.add_argument('-s', '--start-id', type=int,
                        required=False, default=1,
                        dest='start_id',
                        help='Starting allele identifier.')

    args = parser.parse_args()

    return args


if __name__ == '__main__':

    args = parse_arguments()
    main(**vars(args))
