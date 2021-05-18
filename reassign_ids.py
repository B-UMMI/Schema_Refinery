#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 18 17:01:04 2021

@author: rfm
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

    parser.add_argument('-i', type=str, required=True,
                        dest='input_file',
                        help='')

    parser.add_argument('-o', type=str, required=True,
                        dest='output_file',
                        help='')

    parser.add_argument('--s', type=int, required=False,
                        default=1, dest='start_id',
                        help='')

    args = parser.parse_args()

    return args


if __name__ == '__main__':

    args = parse_arguments()
    main(**vars(args))
