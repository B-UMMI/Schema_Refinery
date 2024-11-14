#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import sys
import subprocess
import pandas as pd
from typing import Any, Dict, List

try:
    from DownloadAssemblies import ena661k_assembly_fetcher
    from DownloadAssemblies import ncbi_datasets_summary
    from DownloadAssemblies import ncbi_linked_ids
    from DownloadAssemblies import fetch_metadata
except ModuleNotFoundError:
    from SchemaRefinery.DownloadAssemblies import ena661k_assembly_fetcher
    from SchemaRefinery.DownloadAssemblies import ncbi_datasets_summary
    from SchemaRefinery.DownloadAssemblies import ncbi_linked_ids
    from SchemaRefinery.DownloadAssemblies import fetch_metadata


def find_local_conda_env() -> str:
    """
    Fetches the current conda environment path by using 'conda info'.

    Returns
    -------
    str
        Path to the folder inside the schema refinery environment. 
        May not exist if this module was not used to download assemblies from ENA661K.
    """
    # Run 'conda info' to get information about the conda environment
    conda_path: subprocess.Popen = subprocess.Popen(['conda', 'info'], stdout=subprocess.PIPE)
    stdout, stderr = conda_path.communicate()

    # Decode the output and split it into lines
    conda_info: List[str] = stdout.decode("utf-8").splitlines()

    # Parse the conda info into a dictionary
    info_dict: Dict[str, List[str]] = {}
    for info in conda_info:
        if info == '':
            continue
        key_value: List[str] = info.strip().split(':')
        if len(key_value) == 2:
            info_dict[key_value[0].strip()] = [key_value[1].strip()]
            out_key: str = key_value[0].strip()
        else:
            info_dict[out_key].append(key_value[0].strip())

    # Get the path to the active environment location
    sr_path: str = info_dict['active env location'][0]

    # Return the path to the 'ena661k_files' folder inside the environment
    return os.path.join(sr_path, 'ena661k_files')

def main(args: Any) -> None:
    """
    Main function to handle downloading assemblies based on provided arguments.

    Parameters
    ----------
    args : Any
        Command-line arguments passed to the script.
    """
    # Check for mutually exclusive options
    if args.input_table is not None and args.taxon is not None:
        sys.exit("\nError: Downloading from input table or downloading by taxon are mutually exclusive.")

    # Ensure that either input table or taxon name is provided
    if args.input_table is None and args.taxon is None:
        sys.exit("\nError: Must provide an input table or a taxon name.")

    # Ensure that ENA661K is not used with an input table
    if args.input_table is not None and 'ENA661K' in args.database:
        sys.exit("\nError: Only assemblies from NCBI can be fetched from an input file. ENA661K was parsed.")

    # Handle filtering criteria
    if args.filtering_criteria:
        criteria: Dict[str, Any] = args.filtering_criteria
        if criteria['assembly_source'] == ['GenBank'] and (criteria['verify_status'] is True or criteria['verify_status'] is None):
            sys.exit("\nError: Assembly status can only be verified for assemblies obtained from RefSeq (Set to False, Default(None) = True)")
    else:
        print("\nNo filtering criteria provided.")
        criteria = None

    # Create output directory if it does not exist
    if not os.path.isdir(args.output_directory):
        os.mkdir(args.output_directory)

    print(f"\nFetching assemblies from {args.database} datasets.")
    if 'NCBI' in args.database:
        # Handle NCBI Genome Assembly and Annotation report
        if args.input_table is not None:
            # Read TSV file containing NCBI report
            with open(args.input_table, 'r') as id_list:
                assembly_ids: List[str] = id_list.read().splitlines()

            if len(assembly_ids) == 0:
                print("No assembly identifiers provided.")

            if criteria is not None:
                metadata: Dict[str, Any] = ncbi_datasets_summary.fetch_metadata(args.input_table, None, criteria, args.api_key)
                if metadata['total_count'] == 0:
                    os.sys.exit("\nNo assemblies that satisfy the selected criteria were found.")
        else:
            # Fetch from taxon identifier
            metadata = ncbi_datasets_summary.fetch_metadata(None, args.taxon, criteria, args.api_key)
            if metadata['total_count'] == 0:
                os.sys.exit("\nNo assemblies that satisfy the selected criteria were found.")
            assembly_ids = [sample['accession'] for sample in metadata['reports']]

        total_ids: int = len(assembly_ids)
        failed: List[str] = []
        passed: List[str] = []

        if criteria is not None:
            # Validate assemblies
            if metadata['total_count'] > 0:
                for sample in metadata['reports']:
                    current_accession: str = sample['accession']
                    valid: bool = ncbi_datasets_summary.verify_assembly(sample,
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

        metadata_directory: str = os.path.join(args.output_directory, 'metadata_ncbi')
        if not os.path.exists(metadata_directory):
            os.mkdir(metadata_directory)

        # Save IDs to download
        valid_ids_file: str = os.path.join(metadata_directory, "assemblies_ids_to_download.tsv")
        with open(valid_ids_file, 'w', encoding='utf-8') as ids_to_txt:
            ids_to_txt.write("\n".join(assembly_ids) + '\n')

        # Save IDs that failed criteria
        failed_ids_file: str = os.path.join(metadata_directory, "id_failed_criteria.tsv")
        with open(failed_ids_file, 'w', encoding='utf-8') as ids_to_txt:
            ids_to_txt.write("\n".join(failed) + '\n')

        # If any assembly passed filtering criteria
        if len(assembly_ids) == 0:
            print("\nNo assemblies meet the desired filtering criteria.")
            sys.exit("\nAssemblies that failed are in the following TSV file: {}".format(failed_ids_file))

        # Download assemblies
        if args.download:
            # Build initial arguments for the subprocess run of datasets
            arguments: List[str] = ['datasets', 'download', 'genome', 'accession', '--inputfile', valid_ids_file]

            if args.api_key is not None:
                arguments.extend(['--api-key', args.api_key])

            if criteria is not None and criteria['file_to_include'] is not None:
                arguments.extend(['--include', ','.join(criteria['file_to_include'])])
            else:
                arguments.extend(['--include', 'genome'])

            assemblies_zip: str = os.path.join(args.output_directory, 'assemblies.zip')
            arguments.extend(['--filename', assemblies_zip])
            print("\nDownloading assemblies...")
            subprocess.run(arguments, check=False)
        else:
            print("\nThe list of identifiers for the assemblies that passed the filtering criteria was saved to: {}".format(valid_ids_file))

    # Download from ENA661K
    if 'ENA661K' in args.database:
        # Path for ena661k files
        try:
            # Find local conda path
            sr_path: str = find_local_conda_env()
        except:
            # If not using conda, use output directory instead
            sr_path = os.path.join(args.output_directory, 'ena661k_files')

        metadata_directory = ena661k_assembly_fetcher.main(sr_path,
                                                           args.taxon,
                                                           args.output_directory,
                                                           args.download,
                                                           criteria,
                                                           args.retry,
                                                           args.threads)

    if args.fetch_metadata:
        linked_ids_file: str = os.path.join(metadata_directory, 'id_matches.tsv')
        ids_file: str = os.path.join(metadata_directory, 'assemblies_ids_to_download.tsv')

        print("\nFetching RefSeq, Genbank and SRA IDs linked to the BioSample ID...")
        ncbi_linked_ids.main(ids_file,
                             linked_ids_file,
                             args.email,
                             args.threads,
                             args.retry,
                             args.api_key)

        # Fetch BioSample identifiers
        biosamples: List[str] = pd.read_csv(linked_ids_file, delimiter='\t')['BioSample'].values.tolist()
        # Exclude samples without BioSample identifier
        biosamples = [i for i in biosamples if isinstance(i, str)]
        # Save BioSample identifiers to file
        biosample_file: str = os.path.join(metadata_directory, 'biosamples.tsv')
        with open(biosample_file, 'w+', encoding='utf-8') as ids:
            ids.write("\n".join(biosamples) + '\n')

        print("\nFetching metadata associated to the BioSample ID...")
        fetch_metadata.main(biosample_file,
                            metadata_directory,
                            args.email,
                            args.threads,
                            args.api_key,
                            args.retry)
