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
from typing import Any, Dict, List, Optional


def verify_assembly(metadata_assembly: Dict[str, Any], size_threshold: Optional[float], max_contig_number: Optional[int],
                    genome_size: Optional[int], verify_status: Optional[bool]) -> bool:
    """
    This function verifies assemblies by certain input criteria.

    Parameters
    ----------
    metadata_assembly : Dict[str, Any]
        JSON object (dict) for a single assembly.
    size_threshold : Optional[float]
        (0 >= x >= 1).
    max_contig_number: Optional[int]
        (>0).
    genome_size: Optional[int]
        (>0).
    verify_status: Optional[bool]

    Returns
    -------
    bool
        Boolean value indicating if the assembly passed or failed the criteria.
    """
    # Extract assembly statistics and information
    assembly_stats: Dict[str, Any] = metadata_assembly['assembly_stats']
    assembly_info: Dict[str, Any] = metadata_assembly['assembly_info']

    # Check genome size and size threshold
    if genome_size is not None and size_threshold is not None:
        bot_limit: float = genome_size - (genome_size * size_threshold)
        top_limit: float = genome_size + (genome_size * size_threshold)

        if int(assembly_stats['total_sequence_length']) >= top_limit:
            return False

        if int(assembly_stats['total_sequence_length']) <= bot_limit:
            return False

    # Check maximum contig number
    if max_contig_number is not None:
        if assembly_stats['number_of_contigs'] > max_contig_number:
            return False

    # Check assembly status
    if verify_status is True or verify_status is None:
        if assembly_info['assembly_status'] == 'suppressed':
            return False

    return True


def fetch_metadata(id_list_path: Optional[str], taxon: Optional[str], criteria: Optional[Dict[str, Any]], api_key: Optional[str]) -> Dict[str, Any]:
    """
    This function based on an input id fetches JSON object (dict) for all assemblies.

    Parameters
    ----------
    id_list_path : Optional[str]
        Path to the file containing a list of IDs starting with GCF_ or GCA_.
    taxon : Optional[str]
        Contains desired taxon name.
    criteria: Optional[Dict[str, Any]]
        Contains filtering criteria.
    api_key: Optional[str]
        Key to NCBI API.

    Returns
    -------
    Dict[str, Any]
        JSON object that contains metadata.
    """
    arguments: List[str] = []

    # Add arguments based on id_list_path or taxon
    if id_list_path is not None:
        arguments = ['datasets', 'summary', 'genome', 'accession', '--inputfile', id_list_path]
    elif taxon is not None:
        arguments = ['datasets', 'summary', 'genome', 'taxon', taxon]

    # Add other chosen parameters
    if api_key is not None:
        arguments.extend(['--api-key', api_key])

    if criteria is not None:
        if criteria['assembly_level'] is not None:
            arguments.extend(['--assembly-level', ','.join(criteria['assembly_level'])])
        if criteria['reference'] is True:
            arguments.extend(['--reference'])
        if criteria['exclude_atypical'] is True:
            arguments.extend(['--exclude-atypical'])
        # Filter by chosen assembly source
        if criteria['assembly_source'] is not None:
            arguments.extend(['--assembly-source', ','.join(criteria['assembly_source'])])

    # Run the subprocess to fetch metadata
    metadata_process: subprocess.CompletedProcess = subprocess.run(arguments, stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=False)

    # Parse the JSON output
    metadata: Dict[str, Any] = json.loads(metadata_process.stdout)

    return metadata
