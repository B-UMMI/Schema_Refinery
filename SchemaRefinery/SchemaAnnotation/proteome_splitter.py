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
import gzip
import pickle

from Bio import SeqIO


def proteome_splitter(proteomes_directory: str, output_directory: str):

    # List proteomes in proteome directory
    proteomes = [os.path.join(proteomes_directory, f)
                 for f in os.listdir(proteomes_directory)]

    # Divide into Swiss-Prot and TrEMBL records
    swiss_records = {}
    trembl_records = {}
    descriptions = {}
    for file in proteomes:
        with gzip.open(file, 'rt') as gzfasta:
            for rec in SeqIO.parse(gzfasta, 'fasta'):
                recid = rec.id
                prot = str(rec.seq)
                desc = rec.description
                if recid.startswith('tr'):
                    trembl_records[recid] = prot
                elif recid.startswith('sp'):
                    swiss_records[recid] = prot
                descriptions[recid] = desc

    # Save TrEMBL records to FASTA file
    tr_filename = 'trembl_prots.fasta'
    tr_file_path = os.path.join(output_directory, tr_filename)
    with open(tr_file_path, 'w') as tout:
        records = ['>{0}\n{1}'.format(k, v) for k, v in trembl_records.items()]
        rectext = '\n'.join(records)
        tout.write(rectext+'\n')

    # Save Swiss-Prot records to FASTA file
    sp_filename = 'swiss_prots.fasta'
    sp_file_path = os.path.join(output_directory, sp_filename)
    with open(sp_file_path, 'w') as tout:
        records = ['>{0}\n{1}'.format(k, v) for k, v in swiss_records.items()]
        rectext = '\n'.join(records)
        tout.write(rectext+'\n')

    # Save descriptions with Pickle
    descriptions_filename = 'prots_descriptions'
    descriptions_file_path = os.path.join(output_directory, descriptions_filename)
    with open(descriptions_file_path, 'wb') as dout:
        pickle.dump(descriptions, dout)

    return tr_file_path, sp_file_path, descriptions_file_path
