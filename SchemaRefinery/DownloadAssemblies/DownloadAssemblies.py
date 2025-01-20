#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import subprocess
import pandas as pd
from typing import Any, Dict, List

try:
    from DownloadAssemblies import ena661k_assembly_fetcher
    from DownloadAssemblies import ncbi_datasets
    from DownloadAssemblies import ncbi_linked_ids
    from DownloadAssemblies import fetch_metadata
    from utils import (file_functions as ff,
                       print_functions as pf,
                       logger_functions as logf,
                       globals as gb)
except ModuleNotFoundError:
    from SchemaRefinery.DownloadAssemblies import ena661k_assembly_fetcher
    from SchemaRefinery.DownloadAssemblies import ncbi_datasets
    from SchemaRefinery.DownloadAssemblies import ncbi_linked_ids
    from SchemaRefinery.DownloadAssemblies import fetch_metadata
    from SchemaRefinery.utils import (file_functions as ff,
                                      print_functions as pf,
                                      logger_functions as logf,
                                      globals as gb)


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

def remove_failed_ids(df: pd.DataFrame, failed_ids: List[str]) -> pd.DataFrame:
    """
    Remove rows from a DataFrame that contain failed IDs.

    Parameters
    ----------
    df : pd.DataFrame
        The DataFrame from which to remove rows.
    failed_ids : list of str
        The list of IDs that failed and need to be removed from the DataFrame.

    Returns
    -------
    pd.DataFrame
        A new DataFrame with the rows containing failed IDs removed.

    Examples
    --------
    >>> df = pd.DataFrame({0: ['ID1', 'ID2', 'ID3', 'ID4']})
    >>> failed_ids = ['ID2', 'ID4']
    >>> remove_failed_ids(df, failed_ids)
         0
    0  ID1
    2  ID3
    """
    # Filter the DataFrame to exclude rows where the first column contains any of the failed IDs
    return df[~df[0].isin(failed_ids)]


def main(args: Any) -> None:
    """
    Main function to handle downloading assemblies based on provided arguments.

    Parameters
    ----------
    args : Any
        Command-line arguments passed to the script.
    """
    # Create criteria dictionary
    if args.filtering_criteria:
        criteria: Dict[str, Any] = args.filtering_criteria
    else:
        criteria = None
        pf.print_message("No criteria provided. Fetching all assemblies.", "info")

    # Create output directory if it does not exist
    if not os.path.isdir(args.output_directory):
        os.mkdir(args.output_directory)

    pf.print_message(f"Fetching assemblies from {args.database} datasets.", "info")
    if 'NCBI' in args.database:
        [ncbi_metadata_directory,
         ncbi_valid_ids_file,
         assemblies_zip] = ncbi_datasets.main(args.input_table,
                                                    args.taxon,
                                                    criteria,
                                                    args.output_directory,
                                                    args.download,
                                                    args.api_key)
    else:
        ncbi_metadata_directory = None
        ncbi_valid_ids_file = None
        assemblies_zip = None
    # Download from ENA661K
    if 'ENA661K' in args.database:
        # Path for ena661k files
        try:
            # Find local conda path
            sr_path: str = find_local_conda_env()
        except:
            # If not using conda, use output directory instead
            sr_path = os.path.join(args.output_directory, 'ena661k_files')

        [failed_to_download,
         ena_metadata_directory,
         ena_valid_ids_file,
         assemblies_directory,
         selected_file_ena661k] = ena661k_assembly_fetcher.main(sr_path,
                                                            args.taxon,
                                                            args.output_directory,
                                                            args.download,
                                                            criteria,
                                                            args.retry,
                                                            args.threads)
    else:
        ena_metadata_directory = None
        failed_to_download = []
        ena_valid_ids_file = None
        assemblies_directory = None
        selected_file_ena661k = None

    if args.fetch_metadata:
        all_metadata_directory: str = os.path.join(args.output_directory, 'metadata_all')
        if not os.path.exists(all_metadata_directory):
            ff.create_directory(all_metadata_directory)
        linked_ids_file: str = os.path.join(all_metadata_directory, 'id_matches.tsv')
        ids_file = os.path.join(all_metadata_directory, 'all_ids_fetched.tsv')
        # Merge NCBI and ENA valid IDs
        if ncbi_valid_ids_file and ena_valid_ids_file:
            # Both files are not None, merge them
            ncbi_df = pd.read_csv(ncbi_valid_ids_file, delimiter='\t', header=None)
            ena_df = pd.read_csv(ena_valid_ids_file, delimiter='\t', header=None)
            # Remove failed_to_download IDs
            ena_df = remove_failed_ids(ena_df, failed_to_download)
            # Merge the two dataframes
            merged_df = pd.concat([ncbi_df, ena_df], ignore_index=True)
            # Save to file
            merged_df.to_csv(ids_file, sep='\t', index=False, header=False)
        elif ncbi_valid_ids_file:
            # Only NCBI file is valid
            ncbi_df = pd.read_csv(ncbi_valid_ids_file, delimiter='\t', header=None)
            # Remove failed_to_download IDs
            ncbi_df.to_csv(ids_file, sep='\t', index=False, header=False)
        elif ena_valid_ids_file:
            # Only ENA file is valid
            ena_df = pd.read_csv(ena_valid_ids_file, delimiter='\t', header=None)
            # Remove failed_to_download IDs
            ena_df = remove_failed_ids(ena_df, failed_to_download)
            # Save to file
            ena_df.to_csv(ids_file, sep='\t', index=False, header=False)


        ncbi_valid_ids_file
        pf.print_message("Fetching RefSeq, Genbank and SRA IDs linked to the assembly ID...", "info")

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
        biosample_file: str = os.path.join(all_metadata_directory, 'biosamples_ids.tsv')
        with open(biosample_file, 'w+', encoding='utf-8') as ids:
            ids.write("\n".join(biosamples) + '\n')

        pf.print_message("Fetching metadata associated to the BioSample ID...", "info")

        fetch_metadata.main(biosample_file,
                            all_metadata_directory,
                            args.email,
                            args.threads,
                            args.api_key,
                            args.retry)
    else:
        all_metadata_directory = None
        
    if not args.no_cleanup:
        pf.print_message("Cleaning up temporary files...", "info")
        # Remove temporary files
        ff.cleanup(args.output_directory, [all_metadata_directory,
                                           assemblies_directory,
                                           assemblies_zip,
                                           selected_file_ena661k,
                                           logf.get_log_file_path(gb.LOGGER)])
