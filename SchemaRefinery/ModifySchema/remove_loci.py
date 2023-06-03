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
    
    for loci in loci_to_remove:
        loci_path = os.path.join(schema_path, f"{loci}.fasta")

        if os.path.exists(loci_path):
            # remove locus file from main schema directory
            os.remove(loci_path)
        else:
            print(f"\nERROR: {loci} is not present in the schema, unable to remove.")
