#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
AUTHOR
    Ines Almeida
    github: @IgnesA
DESCRIPTION
    Attempt 1:
    This script accepts a "matches.tsv" table and a "no_match.tsv" 
    table outputed from "match_schemas.py" script. It also accepts one old 
    Loci Annotations table based in an old schema and generates a
    New Annotations Table with the loci ids from the new schema 
    outputted by chewBBACA with the correspondent loci annotations 
    already created.
"""

import argparse
import os
import csv
import copy


def main(input_table_1, input_table_2, input_table_3, output_directory, bsr_threshold):

    processed_matches = []
    dict = {}
    final_table = []

    if not os.path.isdir(output_directory):
        os.mkdir(output_directory)

    

    # open input_table_1 "matches.tsv" outputed from "match_schemas.py" script
    with open(input_table_1, 'r') as table1:
        lines_matches = list(csv.reader(table1, delimiter='\t'))

    # open input_table_2 with loci annotations
    with open(input_table_2, 'r') as table2:
        lines_loci = list(csv.reader(table2, delimiter='\t'))
    
     # open input_table_3 "no_match.tsv" outputed from "match_schemas.py" script
    with open(input_table_3, 'r') as table3:
        lines_no_match = list(csv.reader(table3, delimiter='\t'))

    #processing first table
    for row in lines_matches:
        if ( float(row[2]) >= bsr_threshold ):
            new_list = []
            new_list.append(row[0])

            strs = row[1].split('_')
            new_list.append(strs[0])
            new_list.append(row[2])

            processed_matches.append(new_list)
 

    #processing second table
    for row in lines_loci:
        row[0] = row[0].split('.fasta')[0]
        dict[row[0]] = row

    #creating "new_annotations" table 
    for row in processed_matches:

        entry = dict.get(row[1])

        if entry:
            entry2 = copy.deepcopy(entry)
            entry2[0] = row[0]
            final_table.append(entry2)
        else:
            final_table.append([row[0]])

    #adding the no_match genes to the "new_annotations" table
    for row in lines_no_match:
        final_table.append([row[0]])


    
    with open(output_directory + '/new_annotations.tsv', 'w') as csvfile: 
        # creating a csv writer object 
        csvwriter = csv.writer(csvfile, delimiter='\t') 

        # writing the data rows 
        csvwriter.writerow(lines_loci[0]) #put the collumn names
        csvwriter.writerows(final_table)



    with open(output_directory + '/processed_matches.tsv', 'w') as csvfile: 
        # creating a csv writer object 
        csvwriter = csv.writer(csvfile, delimiter='\t') 

        csvwriter.writerows(processed_matches)

      


def parse_arguments():

    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument('-i1', '--input_table_1', type=str,
                        required=True, dest='input_table_1',
                        help='"matches.tsv" file outputed from "match_schemas.py" script')
    
    parser.add_argument('-i2', '--input_table_2', type=str,
                        required=True, dest='input_table_2',
                        help='TSV file of the old loci annotations table')

    parser.add_argument('-i3', '--input_table_3', type=str,
                        required=True, dest='input_table_3',
                        help='TSV file of the "no_match.txt" file' 
                        'outputed from "match_schemas.py" script')

    parser.add_argument('-bsrt', '--bsr_threshold', type=float,
                        required=False, dest='bsr_threshold', default=0.7,
                        help='Input the minimum bsr threshold acceptable (a number from 0 to 1).' 
                            'Example: 0.6. If this parameter is not given, the default value is 0.7.')

    parser.add_argument('-o', '--output_directory', type=str,
                    required=True, dest='output_directory',
                    help='Path to the directory where generated '
                            'files will be stored')

    args = parser.parse_args()

    return args


if __name__ == '__main__':

    args = parse_arguments()

    main(**vars(args))