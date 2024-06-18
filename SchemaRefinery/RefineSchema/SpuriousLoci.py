#!/usr/bin/env python3
# -*- coding: utf-8 -*-

try:
    from utils import (core_functions as cof)
except ModuleNotFoundError:
    from SchemaRefinery.utils import (core_functions as cof)

def main(schema, output_directory, allelecall_directory, alignment_ratio_threshold, 
        pident_threshold, size_threshold, translation_table, cpu):
    frequency_in_genomes = {}
    loci_ids = [True, True]
    constants = [alignment_ratio_threshold, 
            pident_threshold,
            None,
            None,
            None,
            size_threshold,
            translation_table]

    cof.process_schema(schema,
                       [],
                       output_directory,
                       None,
                       None,
                       frequency_in_genomes,
                       allelecall_directory,
                       None,
                       loci_ids,
                       True,
                       constants,
                       cpu)