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
import hashlib

from Bio import SeqIO

def merge_locus(loci_to_merge,schema_path):
    """Remove the loci contained inside a list
    Parameter
    ---------
    loci_to_remove : list
        List contating id for loci to be removed.

    schema_path : str
        String that contains the new schema path.

    Returns
    -------
    None, operates over OS system folder
    """

    #New loci file name will be the first loci name in the list
    new_allele_id = f">{loci_to_merge[0]}_"
    new_locus_path = os.path.join(schema_path,f"{loci_to_merge[0]}.fasta")
    allele_id = 1
    seq_hash = {}
    failed = False
    for loci in loci_to_merge:
        loci_path = os.path.join(schema_path,f"{loci}.fasta")
        if os.path.exists(loci_path):
            loci_records = [[rec.id, str(rec.seq)] for rec in SeqIO.parse(loci_path, 'fasta')]
            for rec in loci_records:
                hashed_seq = hashlib.sha256(rec[1].encode('utf-8')).hexdigest()
                if hashed_seq in seq_hash.keys():
                    continue
                else:
                    seq_hash[hashed_seq] = [f"{new_allele_id}{allele_id}",
                                            rec[1]]
                    allele_id += 1

            # remove locus file from main schema directory
            os.remove(loci_path)
        else:
            print(f"\nERROR: In the following merge loci group: {' '.join(loci_to_merge)} "
                  f"the following loci : {loci} is not present in the schema. "
                  "Unable to merge with the others")
            failed = True
    
    if failed is False:
        with open(new_locus_path,'w') as outfile:

            for seq in seq_hash.values():
                records = '\n'.join(seq)
                outfile.write(records+'\n')

