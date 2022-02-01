"""
This scripts merges paralagous loci in to groups, and takes output file from from inter_loci_validation.

Inputs:
    -t , path to the inter_loci_validation output file

    -o , path to the output folder

Outputs:
    TSV file containing paralogous groups one per row, each loci is contained in one column.

"""

import os
import argparse
import pandas as pd
import csv
import numpy as np

def main(table_path, output_path):

    table = pd.read_csv(table_path, delimiter="\t")
    table.columns = ['loci','match','a','b','size']

    table = table[table['size'] == 'IN_MODE']

    merged_loci = []
    unique_loci = np.unique([loci.split('_')[0] for loci in pd.unique(table.iloc[0:,0])])
    
    """
    Gets matches to unique loci in the input table e.g:
    x matches y
    x matches z
    x matches k

    return : list containing x,y,z,k and adds it to another list merged_loci containing all groups
    """
    for unique in unique_loci:

        append_list = []
        
        if unique not in merged_loci:
            
            for match in table[table['loci'].str.contains(unique)]['match']:

                if match.split('_')[0] not in append_list:
                    append_list.append(match.split('_')[0])

            append_list.append(unique)  
            append_list.sort()          
            merged_loci.append(append_list)

    """
    As previous for loop returns lists with similiar items there is need to merge similiar lists
    for this reason next for loop finds similiarities in lists and merges them.
    
    e.g :

    [x y z] in [x y z k] deletes first list and maintains the second list

    [x y z] u [x y k] merges to create [x y z k]
    
    """

    for loci in merged_loci:

        for loci2 in merged_loci:

            if loci == loci2:
                continue

            elif set(loci).issubset(set(loci2)):

                if loci in merged_loci:
                    merged_loci.remove(loci)

            elif len(set(loci).intersection(set(loci2)))>0:    

                merged_loci.append(list(set(loci).union(set(loci2))))

                if loci in merged_loci:
                    merged_loci.remove(loci)

                if loci2 in merged_loci:
                    merged_loci.remove(loci2)
    

    for i in range(len(merged_loci)):
        merged_loci[i].sort()

    merged_loci = np.unique(merged_loci)

    with open(os.path.join(output_path,'merged_paralogous_loci.tsv'), 'w', newline='') as f_output:
        tsv_output = csv.writer(f_output, delimiter='\t')

        for i in merged_loci:
            tsv_output.writerow(i)

def parse_arguments():

    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument('-t', type=str, required=True,
                        dest='table_path',
                        help='annotation table')
    
    parser.add_argument('-o', type=str, required=True,
                        dest='output_path',
                        help='output dir')
    

    args = parser.parse_args()

    return args


if __name__ == '__main__':

    args = parse_arguments()

    main(**vars(args))