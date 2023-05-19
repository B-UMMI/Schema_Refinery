#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import sys
import shutil
import csv
from itertools import repeat
from multiprocessing import Pool, cpu_count

try:
    from ModifySchema import merge_loci
except ModuleNotFoundError:
    from SchemaRefinery.ModifySchema import merge_loci

new_schema_path = '/home/ummi/Downloads/test/schema_seed'
input_line = ['merge','a','b']
def multiprocess_table(input_line,new_schema_path):

    if input_line[0].lower() == 'merge':
        merge_loci.merge_loci(input_line[1:],new_schema_path)
    elif input_line[0].lower() == 'remove':
        pass
    elif input_line[0].lower() == 'split':
        pass

def main(args):
    # create output directory
    if os.path.isdir(args.output_directory) is False:
        os.mkdir(args.output_directory)
    
    new_schema_path = os.path.join(args.output_directory,'schema_seed')
    shutil.copytree(args.schema_path, args.new_schema_path)

    with open(args.input_table, 'r', newline='') as file:
        input_table = csv.reader(file, delimiter='\t')

    p = Pool(processes=args.cpu)
    r = p.map_async(multiprocess_table, input_table, repeat(new_schema_path))
    r.wait()



    