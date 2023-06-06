#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import shutil
import concurrent.futures
import subprocess
import pickle

from itertools import repeat

try:
    from ModifySchema import merge_loci
    from ModifySchema import remove_loci
    from ModifySchema import split_loci
except ModuleNotFoundError:
    from SchemaRefinery.ModifySchema import merge_loci
    from SchemaRefinery.ModifySchema import remove_loci
    from SchemaRefinery.ModifySchema import split_loci

def multiprocess_table(input_line,temp_schema_path):
    """Multiprocess
    Parameter
    ---------
    input_line : list
        List containing the command (merge,remove or split) as the first 
        entry followed by loci ids.
    temp_schema_path : str
        String that contains the new schema path.

    Returns
    -------
    None, operates over OS system folder
    """

    if input_line[0].lower() == 'merge':
        merge_loci.merge_locus(input_line[1:],temp_schema_path)
    elif input_line[0].lower() == 'remove':
        remove_loci.remove_locus(input_line[1:],temp_schema_path)
    elif input_line[0].lower() == 'split':
        split_loci.split_locus(input_line[2:],temp_schema_path,input_line[1])

def main(args):

    # create output directory
    if os.path.isdir(args.output_directory) is False:
        print("\nCreating output directory.")
        os.mkdir(args.output_directory)
    
    temp_schema_path = os.path.join(args.output_directory, 'temp_schema_seed')

    terminations = (".fasta",".fna",".fa")
    fasta_files_list = [file for file in os.listdir(args.schema_path) if file.endswith(terminations)]

    if os.path.exists(temp_schema_path):
        print("\ntemp_schema_seed folder already exists...")
        print("Removing and regenerating temp_schema_seed folder.")
        shutil.rmtree(temp_schema_path)
        os.mkdir(temp_schema_path)
        for file in fasta_files_list:
            shutil.copy(file,temp_schema_path)
    else:
        print("\nCopying schema_seed fasta files into temp_schema_seed.")
        os.mkdir(temp_schema_path)
        for file in fasta_files_list:
            from_path = os.path.join(args.schema_path,file)
            to_path = os.path.join(temp_schema_path,file)

            shutil.copy(from_path, to_path)

    print("\nReading input table...")
    input_table = args.input_table

    print("\nExecuting proccesses...")
    with concurrent.futures.ThreadPoolExecutor(max_workers=args.cpu) as executor:
        executor.map(multiprocess_table, input_table, repeat(temp_schema_path))

    print("\nRunning ChewBBACA PrepExternalSchema to regenerate schema_seed.")

    config_file_path = os.path.join(args.schema_path,'.schema_config')

    with open(config_file_path,'rb') as infile:
        config = pickle.load(infile)

    prodigal_training_file = [os.path.join(args.schema_path,file) for file in os.listdir(args.schema_path) if file.endswith('.trn')]

    new_schema_seed_path = os.path.join(args.output_directory,'new_schema_seed')

    if os.path.exists(new_schema_seed_path):
        print("\nnew_schema_seed directory already exists in output directory...")
        print("Removing new_schema_seed directory")
        shutil.rmtree(new_schema_seed_path)

    arguments = ['chewie', 'PrepExternalSchema',
                 '-i', temp_schema_path,
                 '-o', new_schema_seed_path,
                 '--bsr', str(config['bsr'][0]),
                 '--l', str(config['minimum_locus_length'][0]),
                 '--st', str(config['size_threshold'][0]),
                 '--t', str(config['translation_table'][0]),
                 '--ptf', prodigal_training_file[0],
                 '--cpu', str(args.cpu)]

    subprocess.run(arguments)

    print("\nRemoving temp files...")
    os.rmdir(temp_schema_path)
    
    print(f"\nModified Schema available at {new_schema_seed_path}")


    