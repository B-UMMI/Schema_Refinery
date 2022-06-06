#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Author: Inês Almeida
"""

import os
import argparse
import csv


COLUMN_LABEL = 'Paralogous'
OUTPUT_FILE_NAME = "paralogous_groups.tsv"
SEPARATOR = ' | '

def checkAndMakeDirectory(outdir: str):
    if not os.path.isdir(outdir):
        os.mkdir(outdir)

def listTSVFile(file: str):
    with open(file, 'r') as table:
        lines_table = list(csv.reader(table, delimiter='\t'))
    return lines_table

def exportListToTSVFile(filename: str, list: list, first_row: list):
    with open(filename, 'w') as csvfile: 
        # creating a csv writer object 
        csvwriter = csv.writer(csvfile, delimiter='\t') 

        # writing the data rows 
        csvwriter.writerow(first_row) #put the collumn names
        csvwriter.writerows(list)

def get_loci(corr: str):
    return corr.split('&')

def process_input_files(l: list):
    listed_files = [listTSVFile(i) for i in l]
    final = []

    # joining all lines of every file in a single list
    for li in listed_files:
        final += li

    return final

def main(input_files: list, output_directory: str):

    checkAndMakeDirectory(output_directory)

    correspondences_dict = {}
    all_files_merged = process_input_files(input_files)

    for line in all_files_merged:
        loci_list = get_loci(line[1])
        correspondences_found = []
        for loci in loci_list:
            if loci in correspondences_dict.keys():
                correspondences_found.append(loci)
        
        if len(correspondences_found) == 0:
            correspondences = set(loci_list)
            for loci in loci_list:
                correspondences_dict[loci] = correspondences

        else:
            new_set = set()
            for loci in correspondences_found:
                new_set.update(correspondences_dict[loci])
            new_set.update(loci_list)
            for loci in new_set:
                correspondences_dict[loci] = new_set

    unique_correspondences = set()
    for val in correspondences_dict.values():
        s = SEPARATOR.join(val)
        unique_correspondences.add(s)

    exportListToTSVFile(os.path.join(output_directory, OUTPUT_FILE_NAME), [[s] for s in list(unique_correspondences)], [COLUMN_LABEL])    


def parse_arguments():

    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument('-i', type=str, required=True,
                        nargs='+',
                        dest='input_files',
                        help='file to be processed')
    
    parser.add_argument('-o', type=str, required=True,
                        dest='output_directory',
                        help='output directory')
    

    args = parser.parse_args()

    return args


if __name__ == '__main__':

    args = parse_arguments()

    main(**vars(args))