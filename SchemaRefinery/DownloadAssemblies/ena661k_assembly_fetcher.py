#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Purpose
-------
This script enables the download of genome assemblies from the
ENA661k study ("Exploring bacterial diversity via a curated
and searchable snapshot of archived DNA sequences").

Code documentation
------------------
"""

import os
import sys
import csv
import time
import socket
import hashlib
import argparse
import urllib.request
import concurrent.futures
import gzip
import shutil
from itertools import repeat
from typing import List, Tuple, Dict, Any, Union

try:
    from utils import (constants as ct,
                       print_functions as pf)
except ModuleNotFoundError:
    from SchemaRefinery.utils import (constants as ct,
                                      print_functions as pf)


# Increase the field_size_limit to handle large fields in the input file
# "File4_QC_characterisation_661k.txt" passed to this script has fields
# that exceed the default field_size_limit (field: low_coverage_contigs)
maxInt: int = sys.maxsize
while True:
    # This might lead to an OverflowError
    # We need to decrease the value until it is accepted
    try:
        csv.field_size_limit(maxInt)
        break
    except OverflowError:
        maxInt = int(maxInt / 10)

# Set socket timeout for urllib calls to 60 seconds
socket.setdefaulttimeout(60)

def handle_progress(block_num: int, block_size: int, total_size: int) -> None:
    """
    Calculate the progress of the downloads that have been made.

    Parameters
    ----------
    block_num : int
        Number of blocks.
    block_size : int
        Block size in bytes.
    total_size : int
        Size of the file to be downloaded in bytes.

    Returns
    -------
    None
        Prints % of the download in terminal.
    """
    read_data: int = 0
    # Calculate the progress
    # Store a temporary value to store downloaded bytes so that we can
    # add it later to the overall downloaded data
    temp: int = block_num * block_size
    read_data = temp + read_data
    # Calculate the remaining size
    remaining_size: int = total_size - read_data
    if remaining_size <= 0:
        downloaded_percentage: int = 100
        remaining_size = 0
    else:
        downloaded_percentage = int(((total_size - remaining_size) / total_size) * 100)

    if downloaded_percentage == 100:
        pf.print_message(f'Downloaded: {downloaded_percentage}%', "info", end="\n")
    else:
        pf.print_message(" ", None, end="\r")
        pf.print_message(f'Downloaded: {downloaded_percentage}%', "info", end="\r")


def check_download(file: str, file_hash: str, remove: bool = False) -> bool:
    """
    Check the integrity of a downloaded file.

    Parameters
    ----------
    file : str
        Path to the file to check.
    file_hash : str
        Expected md5 hash for the file contents.
    remove : bool, optional
        Whether to remove the file if the hash does not match (default is False).

    Returns
    -------
    bool
        True if the hash determined for the file contents matches
        the expected hash, False otherwise.
    """
    with open(file, 'rb') as infile:
        data: bytes = infile.read()
        md5: str = hashlib.md5(data).hexdigest()

    if md5 != file_hash:
        if remove:
            os.remove(file)
        return False

    return True

def read_table(file_path: str, delimiter: str = '\t') -> List[List[str]]:
    """
    Read a tabular file.

    Parameters
    ----------
    file_path : str
        Path to the tabular file.
    delimiter : str, optional
        Field delimiter (default is '\t').

    Returns
    -------
    List[List[str]]
        List that contains one sublist per line in the input tabular file.
    """
    with open(file_path, 'r', encoding='utf-8') as infile:
        lines: List[List[str]] = list(csv.reader(infile, delimiter=delimiter))

    return lines

def download_ftp_file(data: Tuple[str, str, Union[str, None]], retry: int, verify: bool = True, progress: bool = False) -> bool:
    """
    Download a file from an FTP server.

    Parameters
    ----------
    data : Tuple[str, str, str]
        A tuple containing the FTP path to the file to download, the local path to the file, and the expected md5 hash.
    retry : int
        Maximum number of retries if download fails.
    verify : bool, optional
        Whether to verify the downloaded file's integrity (default is True).
    progress : bool, optional
        Whether to show download progress (default is False).

    Returns
    -------
    bool
        True if the file was successfully downloaded, False otherwise.
    """
    file_url: str = data[0]
    out_file: str = data[1]
    original_hash: str = data[2]
    tries: int = 0
    downloaded: bool = False
    while not downloaded and tries < retry:
        try:
            if not progress:
                urllib.request.urlretrieve(file_url, out_file)
            else:
                urllib.request.urlretrieve(file_url, out_file, handle_progress)
        except Exception:
            time.sleep(1)
            tries += 1
        else:
            # Check if the file was downloaded successfully
            if os.path.isfile(out_file):
                if verify:
                    if check_download(out_file, original_hash, True):
                        downloaded = True
                    else:
                        tries += 1
                else:
                    downloaded = True

    return downloaded

def main(sr_path: str, taxon: str, output_directory: str, ftp_download: bool,
         criteria: Dict[str, Any], retry: int, threads: int) -> tuple[list[str], str, str, str, str]:
    """
    Main function to handle downloading assemblies based on provided arguments.

    Parameters
    ----------
    sr_path : str
        Path for ena661k files.
    taxon : str
        Taxon name to filter the assemblies.
    output_directory : str
        Path to the directory where the output files will be saved.
    ftp_download : bool
        Whether to download files from FTP.
    criteria : Dict[str, Any]
        Filtering criteria for the assemblies.
    retry : int
        Maximum number of retries if download fails.
    threads : int
        Number of threads to use for downloading.

    Returns
    -------
    str
        Path to the metadata directory.
    """
    # FTP paths and local location for files needed to download assemblies
    assembly_ftp_file: str = os.path.join(sr_path, 'assembly_ftp_file.txt')
    assembly_metadata_file: str = os.path.join(sr_path, 'metadata_file.txt')
    local_checklist: str = os.path.join(sr_path, 'checklist.chk')

    # Verify if dir is present in conda env
    if not os.path.exists(sr_path):
        os.mkdir(sr_path)

    # Verify if files are present in the conda dir env
    if os.path.exists(assembly_ftp_file):
        pf.print_message('File with FTP links already exists...', "info")
    else:
        pf.print_message('Downloading ENA661K ftp paths file...', "info")
        download_ftp_file((ct.ASSEMBLY_FTP_PATH, assembly_ftp_file, None), retry, False, True)

    if os.path.exists(assembly_metadata_file):
        pf.print_message('File with ENA661K metadata already exists...', "info")
    else:
        pf.print_message('Downloading ENA661K metadata file...', "info")
        download_ftp_file((ct.ASSEMBLY_METADATA_PATH, assembly_metadata_file + '.gz', None), retry, False, True)

        pf.print_message('Unzipping metadata...', "info")
        with gzip.open(assembly_metadata_file + '.gz', 'rb') as f_in:
            with open(assembly_metadata_file, 'wb') as f_out:
                shutil.copyfileobj(f_in, f_out)

        os.remove(assembly_metadata_file + '.gz')

    if os.path.exists(local_checklist):
        pf.print_message('File with ENA661K checklist already exists...', "info")
    else:
        pf.print_message('Downloading ENA661K checklist.chk...', "info")
        download_ftp_file((ct.FTP_HASH_FILE, local_checklist, None), retry, False, True)

    # Read file with metadata
    pf.print_message("Reading metadata table...", "info")

    metadata_lines: List[List[str]] = read_table(assembly_metadata_file)

    metadata_header: List[str] = metadata_lines[0]

    # Select lines based on taxon name, e.g Brucella, Streptococcus pneumonia
    pf.print_message("Filtering by chosen taxon...", "info")
    taxon_index: int = metadata_header.index('species')
    taxon_lines: List[List[str]] = [line for line in metadata_lines[1:] if all(t in line[taxon_index].split() for t in taxon.split())]

    pf.print_message('\nFound {0} samples for taxon={1}.'.format(len(taxon_lines), taxon), "info")

    if not taxon_lines:
        pf.print_message(f'Did not find matches for {taxon}.', "warning")
        sys.exit(0)

    # Get all ids
    all_sample_ids: List[str] = [line[0] for line in taxon_lines]
    if criteria is not None:
        # Filter based on genome size
        if criteria['genome_size'] is not None and criteria['size_threshold'] is not None:
            bot_limit: float = criteria['genome_size'] - (criteria['genome_size'] * criteria['size_threshold'])
            top_limit: float = criteria['genome_size'] + (criteria['size_threshold'] * criteria['genome_size'])
            size_index: int = metadata_header.index('total_length')
            taxon_lines = [line for line in taxon_lines if int(line[size_index]) >= bot_limit and int(line[size_index]) <= top_limit]

            pf.print_message('{0} with genome size >= {1} and <= {2}.'.format(len(taxon_lines), bot_limit, top_limit), "info")

        # Filter based on taxon abundance
        if criteria['abundance'] is not None:
            abundance_index: int = metadata_header.index('adjust_abundance')
            taxon_lines = [line for line in taxon_lines if float(line[abundance_index]) >= criteria['abundance']]

            pf.print_message('{0} with abundance >= {1}.'.format(len(taxon_lines), criteria['abundance']), "info")

        # Filter based on number of contigs
        if criteria['max_contig_number'] is not None:
            contigs_index: int = metadata_header.index('total_contigs')
            taxon_lines = [line for line in taxon_lines if int(line[contigs_index]) <= criteria['max_contig_number']]

            pf.print_message('{0} with <= {1} contigs.'.format(len(taxon_lines), criteria['max_contig_number']), "info")

        # Filter based on known ST
        if criteria['known_st'] is True:
            st_index: int = metadata_header.index('mlst')
            taxon_lines = [line for line in taxon_lines if line[st_index] != '-']

            pf.print_message('{0} with known ST.'.format(len(taxon_lines)), "info")

        if criteria['ST_list_path'] is not None:
            with open(criteria['ST_list_path'], 'r', encoding='utf-8') as desired_st:
                d_st: List[str] = desired_st.read().splitlines()
                st_index = metadata_header.index('mlst')
                taxon_lines = [line for line in taxon_lines if line[st_index] in d_st]
                pf.print_message('{0} with desired ST'.format(len(taxon_lines)), "info")

        # Filter based on quality level
        if criteria['any_quality'] is False:
            quality_index: int = metadata_header.index('high_quality')
            taxon_lines = [line for line in taxon_lines if line[quality_index] == 'TRUE']

            pf.print_message('{0} with high quality.'.format(len(taxon_lines)), "info")


    # Get sample identifiers
    sample_ids: List[str] = [line[0] for line in taxon_lines]

    # Assemblies that failed filtering criteria
    failed_list: List[str] = [x for x in all_sample_ids if x not in sample_ids]

    ena_metadata_directory: str = os.path.join(output_directory, 'metadata_ena661k')
    if not os.path.exists(ena_metadata_directory):
        os.mkdir(ena_metadata_directory)

    # Write failed and accepted ids to file
    ena_valid_ids_file: str = os.path.join(ena_metadata_directory, "assemblies_ids_to_download.tsv")
    write_type = 'a' if os.path.exists(ena_valid_ids_file) else 'w+'
    with open(ena_valid_ids_file, write_type, encoding='utf-8') as ids_to_tsv:
        ids_to_tsv.write("\n".join(sample_ids) + '\n')

    failed_ids_file: str = os.path.join(ena_metadata_directory, "ids_failed_criteria.tsv")
    with open(failed_ids_file, 'w+', encoding='utf-8') as ids_to_tsv:
        ids_to_tsv.write("\n".join(failed_list) + '\n')

    if len(sample_ids) == 0:
        sys.exit('No assemblies meet the desired filtering criteria.')
    else:
        if criteria is not None:
            pf.print_message('Selected {0} samples/assemblies that meet filtering criteria.'.format(len(sample_ids)), "info")
        else:
            pf.print_message("No filtering criteria were provided. All samples were selected.", "info")

    selected_file_ena661k: str = os.path.join(output_directory, 'assemblies_metadata_ena661k.tsv')
    with open(selected_file_ena661k, 'w', encoding='utf-8') as outfile:
        selected_lines: List[str] = ['\t'.join(line) for line in [metadata_header] + taxon_lines]
        selected_text: str = '\n'.join(selected_lines)
        outfile.write(selected_text + '\n')

    # Putting checksums in dictionary
    hashes_dict: Dict[str, str] = {}
    with open(local_checklist, 'r', encoding='utf-8') as table:
        lines: List[str] = table.readlines()
        for line in lines:
            md5_hash, file_path = line.split('  ')
            file_basename: str = file_path.split('/')[-1].split('.')[0]
            hashes_dict[file_basename] = md5_hash

    if ftp_download:
        # Read table with FTP paths
        ftp_lines: List[List[str]] = read_table(assembly_ftp_file)
        sample_paths: Dict[str, str] = {l[0]: l[1].split('/ebi/ftp')[1] for l in ftp_lines}
        # Get FTP paths only for selected samples
        taxon_paths: Dict[str, str] = {i: sample_paths[i] for i in sample_ids}

        assemblies_directory: str = os.path.join(output_directory, 'ena661k_assemblies')
        if not os.path.exists(assemblies_directory):
            os.mkdir(assemblies_directory)

        # List files in output directory
        local_files: List[str] = os.listdir(assemblies_directory)

        # Create URLs to download
        remote_urls: List[Tuple[str, str, str]] = []
        for sample in sample_ids:
            sample_basename: str = sample_paths[sample].split('/')[-1]
            # Do not download files that have already been downloaded
            if sample_basename not in local_files and sample_basename.split('.gz')[0] not in local_files:
                sample_file: str = os.path.join(output_directory, 'ena661k_assemblies/' + sample_basename)
                sample_url: str = ct.EBI_FTP + sample_paths[sample]
                remote_urls.append((sample_url, sample_file, hashes_dict[sample]))

        if len(remote_urls) < len(sample_ids):
            pf.print_message('{0} assemblies had already been downloaded.'.format(len(sample_ids) - len(remote_urls)), "info")

        pf.print_message('Downloading {0} assemblies...'.format(len(remote_urls)), "info")
        failed: int = 0
        downloaded: int = 0
        with concurrent.futures.ThreadPoolExecutor(max_workers=threads) as executor:
            for res in executor.map(download_ftp_file, remote_urls, repeat(retry)):
                if res:
                    downloaded += 1
                    pf.print_message('Downloaded {0}/{1}'.format(downloaded, len(remote_urls)), "info", end='\r')
                else:
                    failed += 1
        # Print failed downloads
        pf.print_message(f'Failed download for {failed} files.', "info")
        
        failed_to_download = [x for x in sample_ids if x not in [file.split('.')[0] for file in os.listdir(assemblies_directory)]]
        # Write failed downloads to file
        failed_path = os.path.join(ena_metadata_directory, 'failed_to_download.tsv')
        with open(failed_path, 'w', encoding='utf-8') as failed_file:
            failed_file.write('\n'.join(failed_to_download) + '\n')
    else:
        assemblies_directory = None
        failed_to_download = []

    return failed_to_download, ena_metadata_directory, ena_valid_ids_file, assemblies_directory, selected_file_ena661k


def parse_arguments():

    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument('-m', '--metadata-table', type=str,
                        required=True, dest='metadata_table',
                        help='Summarised QC and general characterisation for '
                             'each assembly in the "File4_QC_characterisati'
                             'on_661K" file from the study "Exploring '
                             'bacterial diversity via a curated and '
                             'searchable snapshot of archived DNA '
                             'sequences" (https://doi.org/10.1101/2021.03'
                             '.02.433662).')

    parser.add_argument('-p', '--paths-table', type=str,
                        required=True, dest='paths_table',
                        help='File with sample identifier to FTP path '
                             'mapping (available at http://ftp.ebi.ac.uk/'
                             'pub/databases/ENA2018-bacteria-661k/).')

    parser.add_argument('-s', '--taxon', type=str,
                        required=True, dest='taxon',
                        help='Name of the taxon. Must match one of the '
                             'taxon names in the "taxon" column in the '
                             'metadata table.')

    parser.add_argument('-o', '--output-directory', type=str,
                        required=True, dest='output_directory',
                        help='Path to the output directory.')

    parser.add_argument('--ftp-download', action='store_true',
                        required=False, dest='ftp_download',
                        help='If the assemblies from the selected samples'
                             'should be downloaded.')

    parser.add_argument('-a', '--abundance', type=float,
                        required=False,
                        dest='abundance',
                        help='Minimum taxon abundance. Samples with taxon'
                             ' abundance below this value are not selected.')

    parser.add_argument('-gs', '--genome-size', type=int,
                        required=False,
                        dest='genome_size',
                        help='Expected genome size.')

    parser.add_argument('-st', '--size-threshold', type=float,
                        required=False,
                        dest='size_threshold',
                        help='Genome size can vary in size +/- this value.')

    parser.add_argument('-mc', '--max_contig_number', type=int,
                        required=False,
                        dest='max_contig_number',
                        help='Maximum number of contigs. Assemblies with '
                             'a number of contigs greater than this value '
                             'are not selected.')

    parser.add_argument('--mlst-taxon', type=str,
                        required=False,
                        dest='mlst_taxon',
                        help='The taxon predicted by the MLST tool.')

    parser.add_argument('--known-st', action='store_true',
                        required=False,
                        dest='known_st',
                        help='If the samples must have a known ST.'
                             'Invalid or unkown STs will be "-".')

    parser.add_argument('--any-quality', action='store_true',
                        required=False,
                        dest='any_quality',
                        help='Download all assemblies, even the ones '
                             'that are not high quality.')


    parser.add_argument('-r', '--retry', type=int,
                        required=False, dest='retry',
                        default=7,
                        help='Maximum number of retries when a '
                             'download fails.')

    parser.add_argument('--st', type=str,
                        required=False, dest='st',
                        default=None,
                        help='Desired ST in a txt file,'
                             ' one ST per line')

    args = parser.parse_args()

    return args


if __name__ == "__main__":

    args = parse_arguments()
    main(**vars(args))
