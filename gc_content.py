#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Purpose
-------
Determines global and per sequence GC content for a FASTA file.

Code documentation
------------------
"""


import os
import argparse

from Bio import SeqIO


def main(input_fasta, output_directory):

    # open Fasta file and determine GC content per sequence
    total_bases = {b: 0 for b in ['A', 'T', 'G', 'C']}
    records_gc = {}
    records = SeqIO.parse(input_fasta, 'fasta')
    for rec in records:
        seqid = rec.id
        # get sequence and convert to upper case
        sequence = str(rec.seq).upper()
        A_count = sequence.count('A')
        total_bases['A'] += A_count
        T_count = sequence.count('T')
        total_bases['T'] += T_count
        G_count = sequence.count('G')
        total_bases['G'] += G_count
        C_count = sequence.count('C')
        total_bases['C'] += C_count

        rec_gc = (G_count+C_count) / (G_count+C_count+A_count+T_count)
        records_gc[seqid] = round(rec_gc, 6)

    # calculate gloabal GC content
    global_gc = (total_bases['G']+total_bases['C']) / (total_bases['G']+total_bases['C']+total_bases['A']+total_bases['T'])

    # save results
    input_basename = os.path.basename(input_fasta)
    output_basename = '.'.join(input_basename.split('.')[0:-1])

    detailed_gc = os.path.join(output_directory, output_basename+'_perSeq.tsv')
    with open(detailed_gc, 'w') as outfile:
        detailed_lines = ['{0}\t{1}'.format(k, v)
                          for k, v in records_gc.items()]
        detailed_text = '\n'.join(detailed_lines)
        outfile.write(detailed_text+'\n')

    global_output = os.path.join(output_directory, output_basename+'_global.tsv')
    with open(global_output, 'w') as outfile:
        outfile.write('{0}\t{1}\n'.format(input_basename, round(global_gc, 6)))


def parse_arguments():

    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument('-i', '--input-fasta', type=str,
                        required=True, dest='input_fasta',
                        help='Path to input Fasta file.')

    parser.add_argument('-o', '--output-directory', type=str,
                        required=True, dest='output_directory',
                        help='Path to output directory. Files with '
                             'results are saved to this directory.')

    args = parser.parse_args()

    return args


if __name__ == "__main__":

    args = parse_arguments()
    main(**vars(args))
