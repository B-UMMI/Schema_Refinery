#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Dec 27 13:05:38 2022

AUTHOR

    Mykyta Forofontov
    github: @MForofontov
"""


import json
import subprocess


def verify_assembly(metadata_assembly, size_threshold, max_contig_number,
                    genome_size, verify_status):
    """
    This function verifies assemblies by certain inputa criteria.

    Parameters
    ----------
    metadata_assembly : json object (dict)
        For a single assembly.
    size_threshold : float
        (0 >= x >= 1).
    max_contig_number: int
        (>0).
    genome_size: int
        (>0).
    verify_status: bool

    Returns
    -------
    Boolean value (in order to see if passed or failed)
    """
    assembly_stats = metadata_assembly['assembly_stats']
    assembly_info = metadata_assembly['assembly_info']

    if genome_size is not None and size_threshold is not None:

        bot_limit = genome_size - (genome_size*size_threshold)
        top_limit = genome_size + (genome_size*size_threshold)

        if int(assembly_stats['total_sequence_length']) >= top_limit:

            return False

        if int(assembly_stats['total_sequence_length']) <= bot_limit:

            return False

    if max_contig_number is not None:

        if assembly_stats['number_of_contigs'] > max_contig_number:

            return False

    if verify_status is True or verify_status is None:
        
        if assembly_info['assembly_status'] == 'suppressed':
            return False

    return True


#id_list_path = '/home/rmamede/Desktop/schemaRefinery_PRs/assemblies_ids.txt'
# id_list_path = None
# taxon = taxon
# api_key = None
def fetch_metadata(id_list_path, taxon, criteria, api_key):
    """
    This function based on an input id fetches json object (dict) for all
    assemblies.

    Parameters
    ----------
    input_id : list of str
        starts with GCF_ or GCA_ .
    assembly_level : str
        containing assembly levels separated by ",".
    reference: Bool
    api_key: str

    Returns
    -------
    found_metadata: bool
    failed_list: list
        containing ids that failed initial processing.
    metadata_assembly: json object (dict)
        for all assemblies.
    """
    if id_list_path is not None:
        arguments = ['datasets', 'summary', 'genome', 'accession',
                     '--inputfile', id_list_path]
    elif taxon is not None:
        arguments = ['datasets', 'summary', 'genome', 'taxon', taxon]

    # add other chosen parameters
    if api_key is not None:
        arguments.extend(['--api-key', api_key])

    if criteria is not None:
        if criteria['assembly_level'] is not None:
            arguments.extend(['--assembly-level', ','.join(criteria['assembly_level'])])
        if criteria['reference'] is True:
            arguments.extend(['--reference'])
        if criteria['exclude_atypical'] is True:
            arguments.extend(['--exclude-atypical'])
        # filter by choosen assembly source
        if criteria['assembly_source'] is not None:
            arguments.extend(['--assembly-source', ','.join(criteria['assembly_source'])])

    metadata = subprocess.run(arguments,
                              stdout=subprocess.PIPE,
                              stderr=subprocess.PIPE,
                              check=False)
    print(metadata.stderr)
    metadata = json.loads(metadata.stdout)

    return metadata
