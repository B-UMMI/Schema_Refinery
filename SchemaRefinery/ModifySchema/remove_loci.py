#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Purpose
-------
This script removes loci from a schema.

Code documentation
------------------
"""

import os

def remove_locus(loci_to_remove,schema_path):
    
    for loci in loci_to_remove:
        loci_path = os.path.join(schema_path, f"{loci}.fasta")
        loci_path_short = os.path.join(schema_path, 'short', f"{loci}.fasta")

        if os.path.exists(loci_path) and os.path.exists(loci_path_short):
            # remove locus file from main schema directory
            os.remove(loci_path)
            # remove from short diretory
            os.remove(loci_path_short)
        else:
            print(f"\nERROR: {loci} is not present in the schema, unable to remove.")
