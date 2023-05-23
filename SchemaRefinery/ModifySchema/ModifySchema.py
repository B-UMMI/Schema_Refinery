#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import shutil
import csv
from itertools import repeat
import concurrent.futures

try:
    from ModifySchema import merge_loci
    from ModifySchema import remove_loci
    from ModifySchema import split_loci
except ModuleNotFoundError:
    from SchemaRefinery.ModifySchema import merge_loci
    from SchemaRefinery.ModifySchema import remove_loci
    from SchemaRefinery.ModifySchema import split_loci

def read_tsv(file_path):
    """Read the input TSV and remove the blank columns
    Parameter
    ---------
    file_path : str
        TSV file path.

    Returns
    -------
    data : list
        TSV file converted to list, where each line is a list.
    """
    data = []
    with open(file_path, 'r', newline='') as file:
        reader = csv.reader(file, delimiter='\t')
        for row in reader:
            non_empty_row = [value for value in row if value]
            if non_empty_row:
                data.append(non_empty_row)
    return data

def multiprocess_table(input_line,new_schema_path):
    """Multiprocess
    Parameter
    ---------
    input_line : list
        List containing the command (merge,remove or split) as the first 
        entry followed by loci ids.
    new_schema_path : str
        String that contains the new schema path.

    Returns
    -------
    None, operates over OS system folder
    """

    if input_line[0].lower() == 'merge':
        merge_loci.merge_locus(input_line[1:],new_schema_path)
    elif input_line[0].lower() == 'remove':
        remove_loci.remove_locus(input_line[1:],new_schema_path)
    elif input_line[0].lower() == 'split':
        split_loci.split_locus(input_line[2:],new_schema_path,input_line[1])

def main(args):
    # create output directory
    if os.path.isdir(args.output_directory) is False:
        os.mkdir(args.output_directory)
    
    new_schema_path = os.path.join(args.output_directory,'schema_seed')
    shutil.copytree(args.schema_path, new_schema_path)

    input_table = read_tsv(args.input_table)

    with concurrent.futures.ThreadPoolExecutor(max_workers=args.cpu) as executor:
        executor.map(multiprocess_table, input_table, repeat(new_schema_path))



    