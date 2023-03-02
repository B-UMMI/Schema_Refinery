#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Purpose
-------

This script splits records in a set of UniProt proteomes into
two Fasta files, one with records from Swiss-Prot and another
with records from TrEMBL. It saves the description in the header
of each record into a binary file created with the Pickle module.

Code documentation
------------------
"""

import os
import pickle

from Bio import SeqIO


def proteome_splitter(proteomes_directory:str, output_directory:str):

    print("Starting to split proteomes...")

    # list proteomes in input directory
    proteomes = [os.path.join(proteomes_directory, f)
                 for f in os.listdir(proteomes_directory)
                 if f.endswith('.fasta')]

    # divide into Swiss-Prot and TrEMBL records
    sp = {}
    trembl = {}
    descriptions = {}
    for file in proteomes:
        for rec in SeqIO.parse(file, 'fasta'):
            recid = rec.id
            prot = str(rec.seq)
            desc = rec.description
            if recid.startswith('tr'):
                trembl[recid] = prot
            elif recid.startswith('sp'):
                sp[recid] = prot
            descriptions[recid] = desc

    # save to file
    tr_file = 'trembl_prots.fasta'
    print(f"Saving {tr_file} file...")
    tr_file_path = os.path.join(output_directory, tr_file)
    with open(tr_file_path, 'w') as tout:
        records = ['>{0}\n{1}'.format(k, v) for k, v in trembl.items()]
        rectext = '\n'.join(records)
        tout.write(rectext+'\n')

    sp_file = 'sp_prots.fasta'
    print(f"Saving {sp_file} file...")
    sp_file_path = os.path.join(output_directory, sp_file)
    with open(sp_file_path, 'w') as tout:
        records = ['>{0}\n{1}'.format(k, v) for k, v in sp.items()]
        rectext = '\n'.join(records)
        tout.write(rectext+'\n')

    # save 
    descriptions_file = 'descriptions'
    print(f"Saving {descriptions_file} file...")
    descriptions_file_path = os.path.join(output_directory, descriptions_file)
    with open(descriptions_file_path, 'wb') as dout:
        pickle.dump(descriptions, dout)

    return tr_file_path, sp_file_path, descriptions_file_path
