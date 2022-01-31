#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Purpose
-------
This script merges a group of paralogous loci into a single locus.

Code documentation
------------------
"""


import os
import csv
import shutil
import argparse

from Bio import SeqIO


def main(input_file, schema_dir, output_dir):

    if os.path.isdir(output_dir) is False:
        os.mkdir(output_dir)

    # read groups of loci to merge
    with open(input_file, 'r') as infile:
        loci_groups = list(csv.reader(infile))

    # copy loci to merge to output_dir
    for g in loci_groups:
        for l in g:
            shutil.copy(os.path.join(schema_dir, l+'.fasta'), output_dir)

    merged_loci = {}
    for group in loci_groups:
        loci_ids = group
        # new locus will have the identifier of first locus in list
        new_locus_id = loci_ids[0]
        allele_id = 1
        merged_records = []
        for locus in loci_ids:
            locus_file = os.path.join(output_dir, locus+'.fasta')
            locus_records = [[rec.id, str(rec.seq)] for rec in SeqIO.parse(locus_file, 'fasta')]
            for rec in locus_records:
                if rec[1] not in merged_records:
                    new_id = '>'+rec[0].replace(locus, new_locus_id)
                    new_id = new_id.split('_')
                    new_id[-1] = str(allele_id)
                    new_id = '_'.join(new_id)
                    merged_records.append(new_id)
                    merged_records.append(rec[1])
                    allele_id += 1

        merged_loci[new_locus_id] = merged_records

    # create directory to store merged files
    merged_dir = os.path.join(output_dir, 'merged_loci')
    os.mkdir(merged_dir)

    for k, v in merged_loci.items():
        locus_file = os.path.join(merged_dir, k+'.fasta')
        records = '\n'.join(v)
        with open(locus_file, 'w') as outfile:
            outfile.write(records+'\n')


def parse_arguments():

    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument('-i', type=str, required=True,
                        dest='input_file',
                        help='Path to text file with groups of loci '
                             'to merge. Each line must have a group '
                             'of loci to merge (identifiers separated '
                             'by ",").')

    parser.add_argument('-s', type=str, required=True,
                        dest='schema_dir',
                        help='Path to the schema\'s directory.')

    parser.add_argument('-o', type=str, required=True,
                        dest='output_dir',
                        help='Path to the output directory.')

    args = parser.parse_args()

    return args


if __name__ == '__main__':

    args = parse_arguments()
    main(**vars(args))
