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

def verify_assemblies(metadata,size_threshold, max_contig_number,
                    genome_size, verify_status):
    """
    This function verifies assemblies by certain inputa criteria.

    Parameters
    ----------
    metadata_assembly : json object (dict)
        For all assemblies.
    size_threshold : float
        (0 >= x >= 1).
    max_contig_number: int
        (>0).
    genome_size: int
        (>0).
    verify_status: bool

    Returns
    -------
    list of accepted assemblies metadata.
    """
    
    accepted_assemblies = [meta for meta in metadata["reports"]]
    if genome_size is not None and size_threshold is not None:
        bot_limit = genome_size - (genome_size*size_threshold)
        top_limit = genome_size + (genome_size*size_threshold)

        accepted_assemblies = [meta
                        for meta in accepted_assemblies
                        if int(meta['assembly_stats']['total_sequence_length']) >= bot_limit
                        and int(meta['assembly_stats']['total_sequence_length']) <= top_limit]
        print(f"{len(accepted_assemblies)} with genome size >= {bot_limit} and <= {top_limit}.")

    if max_contig_number is not None:

        accepted_assemblies = [meta
                    for meta in accepted_assemblies
                    if int(meta['assembly_stats']['number_of_contigs']) <= max_contig_number]
        
        print(f"{len(accepted_assemblies)} with <= {max_contig_number} contigs.")
    
    if verify_status is True:

        accepted_assemblies = [meta
                for meta in accepted_assemblies
                if meta['assembly_stats']['assembly_status'] != 'suppressed']
        
        print(f"{len(accepted_assemblies)} with assembly status not suppressed.")
    
    return accepted_assemblies
    

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

    print_out = []
    if id_list_path is not None:
        arguments = ['datasets', 'summary', 'genome', 'accession',
                     '--inputfile', id_list_path]
        print_out.append("Ids")
    elif taxon is not None:
        arguments = ['datasets', 'summary', 'genome', 'taxon', taxon]
        print_out.append(f"Taxon: {taxon}")

    # add other chosen parameters
    if api_key is not None:
        arguments.extend(['--api-key', api_key])

    if criteria is not None:
        if criteria['assembly_level'] is not None:
            arguments.extend(['--assembly-level', ','.join(criteria['assembly_level'])])
            print_out.append(f"with Assembly level: {','.join(criteria['assembly_level'])}")
        if criteria['reference'] is True:
            arguments.extend(['--reference'])
            print_out.append("Reference")
        if criteria['exclude_atypical'] is True or criteria['exclude_atypical'] is None:
            arguments.extend(['--exclude-atypical'])
            print_out.append("Not Atypical")
        # filter by choosen assembly source
        if criteria['assembly_source'] is not None:
            arguments.extend(['--assembly-source', ','.join(criteria['assembly_source'])])
            print_out.append(f"with assembly source: {','.join(criteria['assembly_source'])}")

    metadata = subprocess.run(arguments,
                              stdout=subprocess.PIPE,
                              stderr=subprocess.PIPE,
                              check=False)

    metadata = json.loads(metadata.stdout)

    return metadata,print_out
