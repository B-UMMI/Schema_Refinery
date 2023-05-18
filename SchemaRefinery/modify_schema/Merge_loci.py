#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Purpose
-------
This script merges a group of paralogous loci into a single locus.

Code documentation
------------------
"""


import os
import csv
import shutil
import argparse
import hashlib

from Bio import SeqIO

def merge_loci(loci_to_merge,schema_path):

    #New loci file name will be the first loci name in the list
    new_locus_id = f">{loci_to_merge[0]}_"
    new_locus_path = os.path.join(schema_path,new_locus_path)
    allele_id = 1
    seq_hash = {}

    for loci in loci_to_merge:
        path_to_loci = os.path.join(schema_path,f"{loci}.fasta")

        loci_records = [[rec.id, str(rec.seq)] for rec in SeqIO.parse(path_to_loci, 'fasta')]

        for rec in loci_records:
            hashed_seq = hashlib.sha256(rec[1])

            if hashed_seq in seq_hash.keys():
                continue
            else:
                seq_hash[hashed_seq] = [f"{new_locus_id}{allele_id}",
                                        rec[1]]
                allele_id += 1
        
        os.remove(path_to_loci)
    
    for seq in seq_hash.values():
        records = '\n'.join(seq)

    with open(new_locus_path) as outfile:
        outfile.write(records+'\n')




