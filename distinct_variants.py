#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Purpose
-------
This script determines the distinct variants for a locus based
on the classifications in a matrix of allelic profiles.

Code Documentation
------------------
"""


import os
import argparse

import pandas as pd
from Bio import SeqIO


def main(input_matrix, locus, schema_directory, output_directory):

    matrix_dataframe = pd.read_csv(input_matrix, delimiter='\t')

    # get column with locus classifications
    try:
        locus_classifications = list(matrix_dataframe[locus])
    except Exception as e:
        # column header includes ".fasta" extension
        locus_classifications = list(matrix_dataframe[locus+'.fasta'])

    # get distinct classifications
    distinct_classifications = list(set(locus_classifications))

    # read Fasta file and get sequences that match distinct classifications
    locus_file = os.path.join(schema_directory, locus+'.fasta')

    records = [(rec.id, str(rec.seq)) for rec in SeqIO.parse(locus_file, 'fasta')]

    selected_records = [rec for rec in records
                        if int(rec[0].split('_')[-1]) in distinct_classifications]

    # write file with selected records
    output_file = os.path.join(output_directory, locus+'_distinct.fasta')
    with open(output_file, 'w') as outfile:
        records_lines = ['>{0}\n{1}'.format(*rec) for rec in selected_records]
        records_text = '\n'.join(records_lines)
        outfile.write(records_text+'\n')
    

def parse_arguments():

    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument('-i', type=str, required=True,
                        dest='input_matrix',
                        help='Path to a TSV file that contains a matrix of '
                             'allelic profiles.')

    parser.add_argument('-l', type=str, required=True,
                        dest='locus',
                        help='Locus identifier. Script will determine the '
                             'distinct variants for the locus based on the '
                             'classifications in the input matrix.')

    parser.add_argument('', type=str, required=True,
                        dest='',
                        help='')

    parser.add_argument('-o', type=str, required=True,
                        dest='output_directory',
                        help='Path to the output directory to which the '
                             'Fasta file with the distinct variants will '
                             'be written.')

    args = parser.parse_args()

    return args


if __name__ == '__main__':

    args = parse_arguments()
    main(**vars(args))
