#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Purpose
-------

This script handles everything related with
the creation of TrEMBL and Swiss-Prot annotations.

Code documentation
------------------
"""

from proteome_fetcher import proteomeFetcher
from proteome_splitter import proteomeSplitter
from proteome_matcher import proteomeMatcher

def proteomeAnnotations(input_table:str, proteomes_directory:str, schema_directory:str, threads:int, retry:int, cpu_cores:int, output_directory:str):
    try:
        # proteomes are downloaded to proteomes_directory
        proteomeFetcher(input_table, proteomes_directory, threads, retry)

        tr_file, sp_file, descriptions_file = proteomeSplitter(proteomes_directory, output_directory)

        proteomeMatcher(schema_directory, [tr_file, sp_file, descriptions_file], output_directory, cpu_cores)
    except:
        print("Something went wrong")
