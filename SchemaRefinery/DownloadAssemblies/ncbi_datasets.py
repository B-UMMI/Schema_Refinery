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
import os
import sys
from typing import Any, Dict, List, Optional, Tuple


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

def main(input_table: Optional[str], taxon: Optional[str], criteria: Optional[Dict[str, Any]],
         output_directory: str, download: bool, api_key: Optional[str]) -> Tuple[str, str, str]:
    # Handle NCBI Genome Assembly and Annotation report
    if input_table is not None:
        # Read TSV file containing NCBI report
        with open(input_table, 'r') as id_list:
            assembly_ids: List[str] = id_list.read().splitlines()

        if len(assembly_ids) == 0:
            print("No assembly identifiers provided.")

        if criteria is not None:
            metadata: Dict[str, Any] = fetch_metadata(input_table, None, criteria, api_key)
            if metadata['total_count'] == 0:
                sys.exit("\nNo assemblies that satisfy the selected criteria were found.")
    else:
        # Fetch from taxon identifier
        metadata = fetch_metadata(None, taxon, criteria, api_key)
        if metadata['total_count'] == 0:
            print("\nNo assemblies that satisfy the selected criteria were found for NCBI.")
            continue_run = False
        else:
            continue_run = True
        assembly_ids = [sample['accession'] for sample in metadata['reports']]
    
    if continue_run:
        total_ids: int = len(assembly_ids)
        failed: List[str] = []
        passed: List[str] = []

        if criteria is not None:
            # Validate assemblies
            if metadata['total_count'] > 0:
                for sample in metadata['reports']:
                    current_accession: str = sample['accession']
                    valid: bool = verify_assembly(sample,
                                                criteria['size_threshold'],
                                                criteria['max_contig_number'],
                                                criteria['genome_size'],
                                                criteria['verify_status'])
                    if valid:
                        passed.append(current_accession)
                    else:
                        failed.append(current_accession)
                assembly_ids = passed

            print(f"\n{len(assembly_ids)} passed filtering criteria.")

        ncbi_metadata_directory: str = os.path.join(output_directory, 'metadata_ncbi')
        if not os.path.exists(ncbi_metadata_directory):
            os.mkdir(ncbi_metadata_directory)

        # Save IDs to download
        ncbi_valid_ids_file: str = os.path.join(ncbi_metadata_directory, "assemblies_ids_to_download.tsv")
        with open(ncbi_valid_ids_file, 'w', encoding='utf-8') as ids_to_txt:
            ids_to_txt.write("\n".join(assembly_ids) + '\n')

        # Save IDs that failed criteria
        failed_ids_file: str = os.path.join(ncbi_metadata_directory, "id_failed_criteria.tsv")
        with open(failed_ids_file, 'w', encoding='utf-8') as ids_to_txt:
            ids_to_txt.write("\n".join(failed) + '\n')

        # If any assembly passed filtering criteria
        if len(assembly_ids) == 0:
            print("\nNo assemblies meet the desired filtering criteria.")
            sys.exit("\nAssemblies that failed are in the following TSV file: {}".format(failed_ids_file))

        # Download assemblies
        if download:
            # Build initial arguments for the subprocess run of datasets
            arguments: List[str] = ['datasets', 'download', 'genome', 'accession', '--inputfile', ncbi_valid_ids_file]

            if api_key is not None:
                arguments.extend(['--api-key', api_key])

            if criteria is not None and criteria['file_to_include'] is not None:
                arguments.extend(['--include', ','.join(criteria['file_to_include'])])
            else:
                arguments.extend(['--include', 'genome'])

            assemblies_zip: str = os.path.join(output_directory, 'assemblies_ncbi.zip')
            arguments.extend(['--filename', assemblies_zip])
            print("\nDownloading assemblies...")
            subprocess.run(arguments, check=False)
        else:
            print("\nThe list of identifiers for the assemblies that passed the filtering criteria was saved to: {}".format(ncbi_valid_ids_file))
    
    return ncbi_metadata_directory, ncbi_valid_ids_file, assemblies_zip