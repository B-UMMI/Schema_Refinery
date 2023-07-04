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

try:
    from DownloadAssemblies import constants as ct
except ModuleNotFoundError:
    from SchemaRefinery.DownloadAssemblies import constants as ct


# increase the field_size_limit
# "File4_QC_characterisation_661k.txt" passed to this script has fields
# that exceed field_size_limit (field: low_coverage_contigs)
maxInt = sys.maxsize
while True:
    # this might lead to a OverflowError
    # we need to decrease the value until it is accepted
    try:
        csv.field_size_limit(maxInt)
        break
    except OverflowError:
        maxInt = int(maxInt/10)

# set socket timeout for urllib calls
socket.setdefaulttimeout(60)


# function passed to urllib.request.urlretrieve to track progress
def handle_progress(block_num, block_size, total_size):
    """Calculate %the progress of the downloads that have been made.
    """
    read_data = 0
    # calculating the progress
    # storing a temporary value to store downloaded bytes so that we can
    # add it later to the overall downloaded data
    temp = block_num * block_size
    read_data = temp + read_data
    # calculating the remaining size
    remaining_size = total_size - read_data
    if remaining_size <= 0:
        downloaded_percentage = 100
        remaining_size = 0
    else:
        downloaded_percentage = int(((total_size-remaining_size) / total_size)*(100))

    if downloaded_percentage == 100:
        print(f'Downloaded: {downloaded_percentage}% ', end="\n")
    else:
        print(" ", end="\r")
        print(f'Downloaded: {downloaded_percentage}% ', end="\r")


def check_download(file: str, file_hash: str, remove=False):
    """
    Check the integrity of a downloaded file.

    Parameters
    ----------
    file : str
        Path to the file to check.
    file_hash : str
        Expected md5 hash for the file contents.

    Returns
    -------
    True if the hash determined for the file contents matches
    the expected hash, False otherwise.
    """
    with open(file, 'rb') as infile:
        data = infile.read()
        md5 = hashlib.md5(data).hexdigest()

    if md5 != file_hash:
        if remove is True:
            os.remove(file)
        return False

    return True


def read_table(file_path, delimiter='\t'):
    """
    Read a tabular file.

    Parameters
    ----------
    file_path : str
        Path to the tabular file.
    delimiter : str
        Field delimiter.

    Returns
    -------
    lines : list
        List that contains one sublist per
        line in the input tabular file.
    """
    with open(file_path, 'r', encoding='utf-8') as infile:
        lines = list(csv.reader(infile, delimiter=delimiter))

    return lines


def download_ftp_file(data, retry, verify=True, progress=False):
    """
    Download a file from a FTP server.

    Parameter
    ---------
    file_url : str
        FTP path to the file to download.
    out_file : str
        Local path to the file that will be downloaded.
    retry : int
        Maximum number of retries if download fails.
    original_hash : str

    Returns
    -------
    downloaded : bool
        True if file was successfully downloaded, False otherwise.
    """
    file_url = data[0]
    out_file = data[1]
    original_hash = data[2]
    tries = 0
    downloaded = False
    while downloaded is False and tries < retry:
        try:
            if progress is False:
                res = urllib.request.urlretrieve(file_url, out_file)
            else:
                res = urllib.request.urlretrieve(file_url, out_file, handle_progress)
        except:
            time.sleep(1)
        tries += 1

        if os.path.isfile(out_file) is True:
            if verify is True:
                if check_download(out_file, original_hash, True):
                    downloaded = True
            else:
                downloaded = True

    return downloaded


def main(sr_path, taxon, output_directory, ftp_download, criteria, retry, threads):
    # FTP paths and local location for files needed to download assemblies
    assembly_ftp_file = os.path.join(sr_path, 'assembly_ftp_file.txt')
    assembly_metadata_file = os.path.join(sr_path, 'metadata_file.txt')
    local_checklist = os.path.join(sr_path, 'checklist.chk')

    # Verify if dir is present in conda env
    if not os.path.exists(sr_path):
        os.mkdir(sr_path)

    # Verify if files are present in the conda dir env~
    if os.path.exists(assembly_ftp_file):
        print('\nFile with FTP links already exists...')
    else:
        print('Downloading ENA661K ftp paths file...')
        download_ftp_file([ct.ASSEMBLY_FTP_PATH, assembly_ftp_file, None],
                          retry, False, True)

    if os.path.exists(assembly_metadata_file):
        print('File with ENA661K metadata already exists...')
    else:
        print('Downloading ENA661K metadata file...')
        download_ftp_file([ct.ASSEMBLY_METADATA_PATH, assembly_metadata_file + '.gz', None],
                          retry, False, True)

        print('Unzipping metadata...')
        with gzip.open(assembly_metadata_file + '.gz', 'rb') as f_in:
            with open(assembly_metadata_file, 'wb') as f_out:
                shutil.copyfileobj(f_in, f_out)

        os.remove(assembly_metadata_file + '.gz')

    if os.path.exists(local_checklist):
        print('File with ENA661K checklist already exists...')
    else:
        print('Downloading ENA661K checklist.chk...')
        download_ftp_file([ct.FTP_HASH_FILE, local_checklist, None],
                          retry, False, True)

    # read file with metadata
    print("\nReading metadata table...")
    metadata_lines = read_table(assembly_metadata_file)

    metadata_header = metadata_lines[0]

    # select lines based on taxon name, e.g Brucella, Streptococcus pneumonia
    print("\nFiltering by chosen taxon...")
    taxon_index = metadata_header.index('species')
    taxon_lines = [line
                     for line in metadata_lines[1:]
                     if all(t in line[taxon_index].split() for t in taxon.split())]

    print('\nFound {0} samples for taxon={1}.'
          ''.format(len(taxon_lines), taxon))

    if len(taxon_lines) == 0:
        print(f'Did not find matches for {taxon}.')
        sys.exit(0)

    # get all ids
    all_sample_ids = [line[0] for line in taxon_lines]
    
    if criteria is not None:
        # filter based on genome size
        if criteria['genome_size'] is not None and criteria['size_threshold'] is not None:
            bot_limit = (criteria['genome_size'] -
                        (criteria['genome_size']*criteria['size_threshold']))
            top_limit = (criteria['genome_size'] +
                        (criteria['genome_size']*criteria['size_threshold']))
            size_index = metadata_header.index('total_length')
            taxon_lines = [line
                            for line in taxon_lines
                            if int(line[size_index]) >= bot_limit
                            and int(line[size_index]) <= top_limit]

            print('{0} with genome size >= {1} and '
                '<= {2}.'.format(len(taxon_lines), bot_limit, top_limit))

        # filter based on taxon abundance
        if criteria['abundance'] is not None:
            abundance_index = metadata_header.index('adjust_abundance')
            taxon_lines = [line
                        for line in taxon_lines
                        if float(line[abundance_index]) >= criteria['abundance']]

            print('{0} with abundance >= {1}.'.format(len(taxon_lines),
                                                    criteria['abundance']))

        # filter based on number of contigs
        if criteria['max_contig_number'] is not None:
            contigs_index = metadata_header.index('total_contigs')
            taxon_lines = [line
                        for line in taxon_lines
                        if int(line[contigs_index]) <= criteria['max_contig_number']]

            print('{0} with <= {1} contigs.'.format(len(taxon_lines),
                                                    criteria['max_contig_number']))

        # filter based on known ST
        if criteria['known_st'] is True:
            st_index = metadata_header.index('mlst')
            taxon_lines = [line
                        for line in taxon_lines
                        if line[st_index] != '-']

            print('{0} with known ST.'.format(len(taxon_lines)))

        if criteria['ST_list_path'] is not None:
            with open(criteria['ST_list_path'], 'r', encoding='utf-8') as desired_st:
                d_st = desired_st.read().splitlines()
                st_index = metadata_header.index('mlst')
                taxon_lines = [line
                            for line in taxon_lines
                            if line[st_index] in d_st]
                print('{0} with desired ST'.format(len(taxon_lines)))

        # filter based on quality level
        if criteria['any_quality'] is False:
            quality_index = metadata_header.index('high_quality')
            taxon_lines = [line
                        for line in taxon_lines
                        if line[quality_index] == 'TRUE']

            print('{0} with high quality.'.format(len(taxon_lines)))

    # get sample identifiers
    sample_ids = [line[0] for line in taxon_lines]

    # Assebmlies that failed filtering criteria
    failed_list = [x
                   for x in all_sample_ids
                   if x not in sample_ids]

    metadata_directory = os.path.join(output_directory, 'metadata_ena661k')
    if os.path.exists(metadata_directory) is False:
        os.mkdir(metadata_directory)

    # write failed and accepted ids to file
    valid_ids_file = os.path.join(metadata_directory,
                                  "assemblies_ids_to_download.tsv")
    with open(valid_ids_file, 'w+', encoding='utf-8') as ids_to_tsv:
        ids_to_tsv.write("\n".join(sample_ids)+'\n')

    failed_ids_file = os.path.join(metadata_directory,
                                   "id_failed_criteria.tsv")
    with open(failed_ids_file, 'w+', encoding='utf-8') as ids_to_tsv:
        ids_to_tsv.write("\n".join(failed_list)+'\n')

    if len(sample_ids) == 0:
        sys.exit('\nNo assemblies meet the desired filtering criterias.')
    else:
        print('Selected {0} samples/assemblies that meet filtering '
              'criteria.'.format(len(sample_ids)))

    selected_file = os.path.join(metadata_directory, 'selected_samples.tsv')
    with open(selected_file, 'w', encoding='utf-8') as outfile:
        selected_lines = ['\t'.join(line)
                          for line in [metadata_header]+taxon_lines]
        selected_text = '\n'.join(selected_lines)
        outfile.write(selected_text+'\n')

    # Putting checksums in dictionary
    hashes_dict = {}
    with open(local_checklist, 'r', encoding='utf-8') as table:
        lines = table.readlines()
        for line in lines:
            md5_hash, file_path = line.split('  ')
            file_basename = file_path.split('/')[-1].split('.')[0]
            hashes_dict[file_basename] = md5_hash

    if ftp_download is True:
        # read table with FTP paths
        ftp_lines = read_table(assembly_ftp_file)
        sample_paths = {l[0]: l[1].split('/ebi/ftp')[1] for l in ftp_lines}
        # get FTP paths only for selected samples
        taxon_paths = {i: sample_paths[i] for i in sample_ids}

        assemblies_directory = os.path.join(output_directory, 'ena661k_assemblies')
        if not os.path.exists(assemblies_directory):
            os.mkdir(assemblies_directory)

        # list files in output directory
        local_files = os.listdir(assemblies_directory)

        # create URLs to download
        remote_urls = []
        for sample in sample_ids:
            sample_basename = sample_paths[sample].split('/')[-1]
            # do not download files that have already been downloaded
            if sample_basename not in local_files and sample_basename.split('.gz')[0] not in local_files:
                sample_file = os.path.join(output_directory, 'ena661k_assemblies/' + sample_basename)
                sample_url = ct.EBI_FTP + sample_paths[sample]
                remote_urls.append([sample_url, sample_file, hashes_dict[sample]])

        if len(remote_urls) < len(sample_ids):
            print('{0} assemblies had already been downloaded.'
                  ''.format(len(sample_ids)-len(remote_urls)))

        print('\nDownloading {0} assemblies...'.format(len(remote_urls)))
        failed = 0
        downloaded = 0
        with concurrent.futures.ThreadPoolExecutor(max_workers=threads) as executor:
            for res in executor.map(download_ftp_file, remote_urls, repeat(retry)):
                if res is True:
                    downloaded += 1
                    print('\r', 'Downloaded {0}/{1}'.format(downloaded,
                                                            len(remote_urls)),
                          end='')
                else:
                    failed += 1
            print(f'\nFailed download for {failed} files.')

    return metadata_directory


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
