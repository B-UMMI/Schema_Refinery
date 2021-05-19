#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 19 14:41:38 2021

@author: rfm
"""


import os
import csv
import argparse
import subprocess

from Bio import SeqIO


def main(schema_directory, loci_table, output_directory, adapt_loci):

    # read table with loci to split
    with open(loci_table, 'r') as infile:
        lines = list(csv.reader(infile, delimiter='\t'))

    locus_output_dirs = []
    for line in lines:
        locus = line[0]
        locus_fasta = '{0}.fasta'.format(locus)
        locus_path = os.path.join(schema_directory, locus_fasta)

        locus_output_directory = '{0}/{1}'.format(output_directory, locus+'_out')
        os.mkdir(locus_output_directory)
        locus_output_dirs.append(locus_output_directory)

        # read locus sequences
        sequences = {rec.id: str(rec.seq)
                     for rec in SeqIO.parse(locus_path, 'fasta')}

        # get length intervals
        length_intervals = line[2].split(',')
        length_intervals_int = [i.split('-')
                                for i in length_intervals]
        length_intervals_int = [(int(i[0]), int(i[1])) for i in length_intervals_int]
        length_seqs = {i: [] for i in length_intervals_int}

        # select sequences for each length interval
        for k in length_seqs:
            min_length = k[0]
            max_length = k[1]
            length_seqs[k] = [recid
                              for recid, seq in sequences.items()
                              if (min_length <= len(seq) <= max_length)]

        # save sequences to new files
        new_ids = line[1].split(',')
        new_ids_map = {k: new_ids[i] for i, k in enumerate(length_seqs.keys())}

        for k, v in new_ids_map.items():
            output_file = '{0}/{1}.fasta'.format(locus_output_directory, v)
            records = {v+'_'+str(i+1): sequences[recid]
                       for i, recid in enumerate(length_seqs[k])}

            # save records
            fasta_records = ['>{0}\n{1}'.format(recid, seq)
                             for recid, seq in records.items()]
            fasta_text = '\n'.join(fasta_records)
            with open(output_file, 'w') as outfile:
                outfile.write(fasta_text+'\n')

    # create output with mapping between old and new sequences!

    # adapt with PrepExternalSchema
    if adapt_loci is True:
        for d in locus_output_dirs:
            adapted_dir = os.path.join(d, 'adapted')
            prep_res = subprocess.Popen(['chewBBACA.py', 'PrepExternalSchema', '-i',
                                         d, '-o', adapted_dir, '--st', 'None'],
                                        stdout=subprocess.PIPE,
                                        stderr=subprocess.PIPE)

            stderr = prep_res.stderr.readlines()

            if len(stderr) > 0:
                print(stderr)


def parse_arguments():

    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument('-s', type=str, required=True,
                        dest='schema_directory',
                        help='')

    parser.add_argument('-t', type=str, required=True,
                        dest='loci_table',
                        help='')

    parser.add_argument('-o', type=str, required=True,
                        dest='output_directory',
                        help='')

    parser.add_argument('--adapt', required=False, action='store_true',
                        dest='adapt_loci',
                        help='')

    args = parser.parse_args()

    return args


if __name__ == '__main__':

    args = parse_arguments()
    main(**vars(args))
