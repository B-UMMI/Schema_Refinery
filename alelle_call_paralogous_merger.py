"""
This script merges paralogous loci from RepeatedLoci.txt in alelle_call results.

Inputs :
    -t , takes various paths to RepeatedLoci.txt in alelle calls, e.g path_1,path_2,path_3,path_4
    -s , path to the schema
    -o , output path

Outputs :
    TSV file containing paralogous groups one per row, each loci is contained in one column.

"""

import os
import argparse
import pandas as pd
import csv
import numpy as np
from collections import Counter
from Bio import SeqIO
import statistics
import math

def mean_dict_calc(locus_path):

    """
    Calculates mean loci size:

    Arguments:
        locus_path : str
            path to fasta file

    Returns :
        locus_mean : int
            mean value for the following loci
    """

    records = SeqIO.parse(locus_path, 'fasta')
    seqs_dict = {rec.id: str(rec.seq.upper()) for rec in records}
    
    # get length of each allele
    allele_lengths = [len(seq) for seq in seqs_dict.values()]

    # determine sequence length mean
    locus_mean = round(statistics.mean(allele_lengths), 1)

    return locus_mean


def merger(table_path):

    """
    Merges items by same PC values:

    Arguments:
        table_path : str
            path to RepeatedLoci.txt

    Returns :
        merged_loci : list
            list of locis merged by PC value
    """

    table = pd.read_csv(table_path, delimiter="\t")

    unique_cp = pd.unique(table['PC'])
    merged_loci = []

    for unique in unique_cp:
        append_list = []
        for id in table.loc[table['PC'] == unique]['gene']:
            append_list.append(id)

        append_list.sort()

        merged_loci.append(append_list)

    merged_loci.sort()
    merged_loci = np.unique(merged_loci)

    return merged_loci

def main(table_path, schema_path, output_path):

    fastas = [f for f in os.listdir(schema_path) if '.fasta' in f]

    fastas = [os.path.join(schema_path, f) for f in fastas]

    mean_dict = {}

    for locus in fastas:
        mean_dict[locus.split('/')[-1]] = mean_dict_calc(locus)


    merged_table_all = []

    """
    Obtain similar results between different paralagous output files and adds them
    to common table.

    if all elements of one selected group are with similiar mean sizes(<=20% size difference)
    they are kept otherwise they are removed.
    """

    for table in table_path.split(','):

        merged_table = merger(table)

        for loci_group in merged_table:

            loci_group_bool = []

            for loci1 in loci_group:

                for loci2 in loci_group:

                    if loci1 == loci2:

                        continue

                    else:
                        loci_group_bool.append(math.isclose(mean_dict[loci1],mean_dict[loci2],rel_tol = 0.2))
            
            if all(loci_group_bool):
                merged_table_all.append(loci_group)

    
    """
    Verifies if two groups have atleast one similiar element, if they have following code merges them
    and chooses only unique items in merged_table_all variable.

    e.g :

    [x y z] in [x y z k] deletes first list and maintains the second list

    [x y z] u [x y k] merges to create [x y z k]

    """

    loop = 0
    loops_todo = Counter([i[0] for i in merged_table_all]).most_common()[0][1]

    while loops_todo != loop:

        loop +=1

        for loci_group in merged_table_all:

            for loci_group_2 in merged_table_all:

                if len(set(loci_group).intersection(set(loci_group_2)))>0:    

                    merged_table_all.append(list(set(loci_group).union(set(loci_group_2))))

                    if loci_group_2 in merged_table_all:

                        merged_table_all.remove(loci_group_2)

        merged_table_all = list(np.unique(merged_table_all))                

    for i in range(len(merged_table_all)):

        merged_table_all[i].sort()

    merged_table_all = np.unique(merged_table_all)


    with open(os.path.join(output_path,'merged_paralog_loci.tsv'), 'w', newline='') as f_output:

        tsv_output = csv.writer(f_output, delimiter='\t')

        for loci in merged_table_all:

            tsv_output.writerow(loci)

def parse_arguments():

    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument('-t', type=str, required=True,
                        dest='table_path',
                        help='allele call parallagous result path separated by ,')

    parser.add_argument('-s', type=str, required=True,
                        dest='schema_path',
                        help='schema_path')
    
    parser.add_argument('-o', type=str, required=True,
                        dest='output_path',
                        help='output dir')
    

    args = parser.parse_args()

    return args


if __name__ == '__main__':

    args = parse_arguments()

    main(**vars(args))

