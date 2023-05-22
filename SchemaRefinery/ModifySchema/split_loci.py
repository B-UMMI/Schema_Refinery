#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os

from Bio import SeqIO


def split_locus(loci_list, schema_path, loci_to_split):
    loci_path = os.path.join(schema_path,f"{loci_to_split}.fasta")
    loci_path_short = os.path.join(schema_path, 'short', f"{loci_to_split}.fasta")

    if os.path.exists(loci_path) and os.path.exists(loci_path_short):
        sequences = {rec.id: str(rec.seq)
                    for rec in SeqIO.parse(locus_path, 'fasta')}

    for new_loci_name, interval in zip(*[iter(loci_list)]*2):
        min_length, max_length = interval.split('-')

        new_loci_path = os.path.join(schema_path,f"{new_loci_name}.fasta")

        new_allele_id = f">{new_loci_name}_"

        new_allele_list = [seq for seq in sequences.items() 
                         if len(seq)>=min_length or len(seq)<=max_length]
        allele_id = 0
        new_allele_dict = {}
        for allele in new_allele_list:
            new_allele_dict[f"{new_allele_id}{allele_id}"] = [f"{new_allele_id}{allele_id}",allele]

        with open(new_loci_path,'w') as outfile:

            for seq in new_allele_dict.values():
                records = '\n'.join(seq)
                outfile.write(records+'\n')
        
    # remove locus file from main schema directory
    os.remove(loci_path)
    # remove from short diretory
    os.remove(loci_path_short)
