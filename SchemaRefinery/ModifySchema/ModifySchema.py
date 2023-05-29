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

    print("Creating directory and copying schema...")
    # create output directory
    if os.path.isdir(args.output_directory) is False:
        os.mkdir(args.output_directory)
    
    new_schema_path = os.path.join(args.output_directory,'schema_seed')

    if os.path.exists(new_schema_path):
        shutil.rmtree(new_schema_path)
        shutil.copytree(args.schema_path, new_schema_path)
    else:
        shutil.copytree(args.schema_path, new_schema_path)

    print("Reading input table...")
    input_table = args.input_table

    print("Executing proccesses...")
    with concurrent.futures.ThreadPoolExecutor(max_workers=args.cpu) as executor:
        executor.map(multiprocess_table, input_table, repeat(new_schema_path))

    print(f"Modified Schema available at {new_schema_path}")


    