#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

@author: rfm
"""


import os
import csv
import argparse

from Bio import SeqIO


def main(schema, locus_id, allelecall_matrix, output_dir):

    if os.path.isdir(output_dir) is False:
        os.mkdir(output_dir)

    # import AlleleCall matrix
    with open(allelecall_matrix, 'r') as infile:
        profiles = list(csv.reader(infile, delimiter='\t'))

    # get locus column
    locus_index = profiles[0].index(locus_id+'.fasta')

    # get locus classifications for all genomes
    locus_classifications = {}
    for line in profiles[1:]:
        allele_id = line[locus_index]
        if 'INF-' in allele_id:
            allele_id = allele_id.replace('INF-', '')

        try:
            allele_id = int(allele_id)
            sample_id = '.'.join(line[0].split('.')[:-1])
            locus_classifications.setdefault(allele_id, []).append(sample_id)
        except ValueError:
            pass

    # determine frequency of each allele
    for k, v in locus_classifications.items():
        v.append(len(v))

    # determine alleles that were not in any genome
    locus_file = os.path.join(schema, locus_id+'.fasta')
    alleles_ids = [int((rec.id).split('_')[-1])
                   for rec in SeqIO.parse(locus_file, 'fasta')]

    for a in alleles_ids:
        if a not in locus_classifications:
            locus_classifications[a] = ['', 0]

    # create output file lines
    output_lines = [[k, v[-1], v[:-1]] for k, v in locus_classifications.items()]
    # sort based on decreased frequency
    output_lines = sorted(output_lines, key=lambda x: x[1], reverse=True)
    output_header = 'allele_id\tfrequency\tsamples'
    output_text = ['{0}\t{1}\t{2}'.format(l[0], l[1], ','.join(l[2]))
                   for l in output_lines]
    output_text = [output_header] + output_text
    output_text = '\n'.join(output_text)

    output_file = os.path.join(output_dir, '{0}_frequencies.tsv'.format(locus_id))
    with open(output_file, 'w') as outfile:
        outfile.write(output_text)


def parse_arguments():

    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument('-s', type=str, required=True,
                        dest='schema',
                        help='')

    parser.add_argument('-l', type=str, required=True,
                        dest='locus_id',
                        help='')

    parser.add_argument('-m', type=str, required=True,
                        dest='allelecall_matrix',
                        help='')

    parser.add_argument('-o', type=str, required=True,
                        dest='output_dir',
                        help='')

    args = parser.parse_args()

    return args


if __name__ == '__main__':

    args = parse_arguments()

    main(**vars(args))
