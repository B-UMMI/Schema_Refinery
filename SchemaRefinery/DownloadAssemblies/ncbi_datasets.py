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
from typing import Any, Dict, List, Optional, Tuple, Union

try:
    from utils import (Types as tp,
                       print_functions as pf)
except ModuleNotFoundError:
    from SchemaRefinery.utils import (Types as tp,
                                      print_functions as pf)

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


def fetch_metadata(id_list_path: Optional[str], taxon: Optional[str], criteria: Optional[Dict[str, Any]],
                   api_key: Optional[str]) -> tp.Metadata:
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
    tp.Metadata
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
    metadata: tp.Metadata = json.loads(metadata_process.stdout)

    return metadata

def main(input_table: Optional[str], taxon: Optional[str], criteria: Optional[Dict[str, Any]],
         output_directory: str, download: bool, api_key: Optional[str]) -> Tuple[str, str, Union[str, None]]:
    # Handle NCBI Genome Assembly and Annotation report
    if input_table is not None:
        # Read TSV file containing NCBI report
        with open(input_table, 'r') as id_list:
            assembly_ids: List[str] = id_list.read().splitlines()

        if len(assembly_ids) == 0:
            pf.print_message("No assembly IDs were found in the input file.", "warning")

        if criteria is not None:
            metadata: tp.Metadata = fetch_metadata(input_table, None, criteria, api_key)
            if metadata['total_count'] == 0:
                pf.print_message("No assemblies that satisfy the selected criteria were found for NCBI.", "warning")
                sys.exit()
    else:
        # Fetch from taxon identifier
        metadata = fetch_metadata(None, taxon, criteria, api_key)
        if metadata['total_count'] == 0:
            pf.print_message("No assemblies that satisfy the selected criteria were found for NCBI.", "warning")
            continue_run = False
        else:
            continue_run = True
        assembly_ids = [sample['accession'] for sample in metadata['reports']]

        assemblies_zip: Union[str, None] = None
    
    if continue_run:
        total_ids: int = len(assembly_ids)
        failed: List[str] = []
        passed: List[str] = []
        passed_metadata: List[Dict[str, Any]] = []
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
                        # Append to passed list
                        passed.append(current_accession)
                        # Append metadata to list
                        passed_metadata.append(sample)
                    else:
                        failed.append(current_accession)
                # Update assembly_ids to only include passed assemblies
                assembly_ids = passed
                # Update metadata to only include passed assemblies
                metadata = {'total_count': len(assembly_ids), 'reports': passed_metadata}

            pf.print_message(f"{len(assembly_ids)} passed filtering criteria.", "info")

        ncbi_metadata_directory: str = os.path.join(output_directory, 'metadata_ncbi')
        if not os.path.exists(ncbi_metadata_directory):
            os.mkdir(ncbi_metadata_directory)

        metadata_file: str = os.path.join(output_directory, 'assemblies_metadata_ncbi.tsv')
        with open(metadata_file, 'w', encoding='utf-8') as metadata_tsv:
            # Write header
            metadata_tsv.write(
                                "ID\t"
                                "assembly_level\t"
                                "biosample_ID\t"
                                "bioproject_ID\t"
                                "taxon\t"
                                "submission_date\t"
                                "contig_l50\t"
                                "contig_n50\t"
                                "gc_count\t"
                                "gc_percent\t"
                                "genome_coverage\t"
                                "number_of_component_sequences\t"
                                "number_of_contigs\t"
                                "number_of_scaffolds\t"
                                "scaffold_l50\t"
                                "scaffold_n50\t"
                                "total_number_of_chromosomes\t"
                                "total_sequence_length\t"
                                "total_ungapped_length\n"
            )
            # Loop through metadata and write to file
            # HINT: you can see all available keys the whole set of assemblies by using itf.identify_dict_structure(metadata[reports])
            for sample in metadata['reports']:
                assembly_info = sample.get('assembly_info', {})
                biosample = assembly_info.get('biosample', {})
                assembly_stats = sample.get('assembly_stats', {})
                data = {
                    "ID": sample.get('accession', 'NA'),
                    "assembly_level": assembly_info.get('assembly_level', 'NA'),
                    "biosample_ID": biosample.get('accession', 'NA'),
                    "bioproject_IDs": assembly_info.get('bioproject_accession', 'NA'),
                    "taxon": sample.get('organism', {}).get('organism_name', 'NA'),
                    "submission_date": biosample.get('submission_date', 'NA'),
                    "contig_l50": assembly_stats.get('contig_l50', 'NA'),
                    "contig_n50": assembly_stats.get('contig_n50', 'NA'),
                    "gc_count": assembly_stats.get('gc_count', 'NA'),
                    "gc_percent": assembly_stats.get('gc_percent', 'NA'),
                    "genome_coverage": assembly_stats.get('genome_coverage', 'NA'),
                    "number_of_component_sequences": assembly_stats.get('number_of_component_sequences', 'NA'),
                    "number_of_contigs": assembly_stats.get('number_of_contigs', 'NA'),
                    "number_of_scaffolds": assembly_stats.get('number_of_scaffolds', 'NA'),
                    "scaffold_l50": assembly_stats.get('scaffold_l50', 'NA'),
                    "scaffold_n50": assembly_stats.get('scaffold_n50', 'NA'),
                    "total_number_of_chromosomes": assembly_stats.get('total_number_of_chromosomes', 'NA'),
                    "total_sequence_length": assembly_stats.get('total_sequence_length', 'NA'),
                    "total_ungapped_length": assembly_stats.get('total_ungapped_length', 'NA')
                }
                # Write data to file
                metadata_tsv.write(
                                f"{data['ID']}\t"
                                f"{data['assembly_level']}\t"
                                f"{data['biosample_ID']}\t"
                                f"{data['bioproject_IDs']}\t"
                                f"{data['taxon']}\t"
                                f"{data['submission_date']}\t"
                                f"{data['contig_l50']}\t"
                                f"{data['contig_n50']}\t"
                                f"{data['gc_count']}\t"
                                f"{data['gc_percent']}\t"
                                f"{data['genome_coverage']}\t"
                                f"{data['number_of_component_sequences']}\t"
                                f"{data['number_of_contigs']}\t"
                                f"{data['number_of_scaffolds']}\t"
                                f"{data['scaffold_l50']}\t"
                                f"{data['scaffold_n50']}\t"
                                f"{data['total_number_of_chromosomes']}\t"
                                f"{data['total_sequence_length']}\t"
                                f"{data['total_ungapped_length']}\n"
                            )

        # Save IDs to download
        ncbi_valid_ids_file: str = os.path.join(ncbi_metadata_directory, "assemblies_ids_to_download.tsv")
        with open(ncbi_valid_ids_file, 'w', encoding='utf-8') as ids_to_txt:
            ids_to_txt.write("\n".join(assembly_ids) + '\n')

        # Save IDs that failed criteria
        failed_ids_file: str = os.path.join(ncbi_metadata_directory, "ids_failed_criteria.tsv")
        with open(failed_ids_file, 'w', encoding='utf-8') as ids_to_txt:
            ids_to_txt.write("\n".join(failed) + '\n')

        # If any assembly passed filtering criteria
        if len(assembly_ids) == 0:
            pf.print_message("No assemblies meet the desired filtering criteria.", "warning")
            pf.print_message("Assemblies that failed are in the following TSV file: {}".format(failed_ids_file), "info")
            sys.exit()

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

            assemblies_zip = os.path.join(output_directory, 'assemblies_ncbi.zip')
            arguments.extend(['--filename', assemblies_zip])
            pf.print_message("Downloading assemblies...", "info")
            subprocess.run(arguments, check=False)
        else:
            pf.print_message("The list of identifiers for the assemblies that passed the filtering criteria was saved to: {}".format(ncbi_valid_ids_file), "info")

    return ncbi_metadata_directory, ncbi_valid_ids_file, assemblies_zip