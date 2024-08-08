#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import sys

try:
    from utils import (core_functions as cof)
except ModuleNotFoundError:
    from SchemaRefinery.utils import (core_functions as cof)
    
def main(schema, output_directory, allelecall_directory, alignment_ratio_threshold_gene_fusions, 
        pident_threshold_gene_fusions, clustering_sim, clustering_cov, genome_presence,
        size_threshold, translation_table, bsr, problematic_proportion, cpu):
    
    temp_paths = [os.path.join(allelecall_directory, "temp"), 
                      os.path.join(allelecall_directory, "unclassified_sequences.fasta"),
                      os.path.join(allelecall_directory, "missing_classes.fasta")]
    # Put all constants in one dict in order to decrease number of variables
    # used around.
    constants = [alignment_ratio_threshold_gene_fusions, 
                pident_threshold_gene_fusions,
                genome_presence,
                clustering_sim,
                clustering_cov,
                size_threshold,
                translation_table,
                bsr,
                problematic_proportion]
    
    if not os.path.exists(temp_paths[0]) or not os.path.exists(temp_paths[1]):
        sys.exit(f"Error: {temp_paths[0]} must exist, make sure that AlleleCall "
                    "was run using --no-cleanup and --output-unclassified flag.")

    unclassified_cds_output = os.path.join(output_directory, "unclassified_cds")
    cof.classify_cds(schema, unclassified_cds_output, allelecall_directory,
                constants, temp_paths, cpu)