#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Purpose
-------
This script accepts a Genome Assembly and Annotation report
table from the NCBI and downloads the genome assemblies of
the samples listed in the table.

Code documentation
------------------
"""


import os
import sys
import csv
import time
import socket
import argparse
import concurrent.futures
from itertools import repeat
from typing import List

try:
    from utils import download_functions as df
except:
    from SchemaRefinery.utils import download_functions as df

# Set socket timeout for urllib calls
socket.setdefaulttimeout(30)

# URL to download assembly_summary_refseq.txt
assembly_summary_refseq: str = 'https://ftp.ncbi.nlm.nih.gov/genomes/ASSEMBLY_REPORTS/assembly_summary_refseq.txt'

def main(input_table: str, output_directory: str, file_extension: str,
         ftp: List[str], threads: int, species: str, retry: int) -> None:
    """
    Main function to handle downloading assemblies based on provided arguments.

    Parameters
    ----------
    input_table : str
        Path to the input table downloaded from NCBI.
    output_directory : str
        Path to the directory where the output files will be saved.
    file_extension : str
        File extension for the files to be downloaded.
    ftp : List[str]
        List of FTP sources to download from (e.g., ['refseq', 'genbank']).
    threads : int
        Number of threads to use for downloading.
    species : str
        Species name to filter the assemblies.
    retry : int
        Maximum number of retries if download fails.

    Returns
    -------
    None
        The function writes the output files to the specified directory.
    """

    species = species.lower()

    # Create output directory if it does not exist
    if not os.path.isdir(output_directory):
        os.mkdir(output_directory)

    # Open table downloaded from NCBI
    with open(input_table, 'r') as table:
        lines: List[List[str]] = list(csv.reader(table, delimiter=','))

    # Get URLs for samples that have RefSeq FTP path
    refseq_urls: List[str] = []
    refseq_assemblies_ids: List[str] = []
    if 'refseq' in ftp:
        refseq_urls = [line[15] for line in lines[1:] if line[15].strip() != '']
        refseq_assemblies_ids = [url.split('/')[-1] for url in refseq_urls]
        print('Found URLs for {0} strains in RefSeq.'.format(len(refseq_urls)))

    genbank_urls: List[str] = []
    genbank_assemblies_ids: List[str] = []
    if 'genbank' in ftp:
        # Get GenBank URLs, excluding samples that are also in RefSeq
        genbank_urls = [line[14] for line in lines[1:] if line[14].strip() != '' and line[15] not in refseq_urls]
        genbank_assemblies_ids = [url.split('/')[-1] for url in genbank_urls]
        print('Found URLs for {0} strains in GenBank.'.format(len(genbank_urls)))

    print('Downloading assembly_summary_refseq...', end='')
    assembly_summary_refseq_local: str = os.path.join(output_directory, 'assembly_summary_refseq.txt')
    df.download_file(assembly_summary_refseq, assembly_summary_refseq_local, retry)
    print('done.')

    # Open assembly_summary_refseq table
    with open(assembly_summary_refseq_local, 'r') as table1:
        assembly_summary_lines: List[List[str]] = list(csv.reader(table1, delimiter='\t'))

    os.remove(assembly_summary_refseq_local)

    # Process assembly_summary_refseq table, filtering for the species given as argument
    assembly_latest: List[str] = []
    for row in assembly_summary_lines[1:]:
        if species in row[7].lower():
            assembly_latest.append(row[0])

    # Remove suppressed genomes from refseq_assemblies_ids and refseq_urls
    suppressed: List[int] = []
    for i in refseq_assemblies_ids:
        processed_id: str = '_'.join(i.split('_')[0:2])
        if processed_id not in assembly_latest:
            suppressed.append(refseq_assemblies_ids.index(i))

    print('Found {0} ids suppressed from RefSeq.'.format(len(suppressed)))

    # Remove suppressed from lists
    refseq_assemblies_ids = [j for i, j in enumerate(refseq_assemblies_ids) if i not in suppressed]
    refseq_urls = [j for i, j in enumerate(refseq_urls) if i not in suppressed]

    urls: List[str] = refseq_urls + genbank_urls
    assemblies_ids: List[str] = refseq_assemblies_ids + genbank_assemblies_ids

    # Construct FTP URLs
    ftp_urls: List[str] = []
    for i, url in enumerate(urls):
        ftp_url: str = '{0}/{1}_{2}'.format(url, assemblies_ids[i], file_extension)
        ftp_urls.append(ftp_url)

    files_number: int = len(ftp_urls)
    if files_number == 0:
        sys.exit('No valid ftp links after scanning table.')

    assemblies_ids = ['{0}/{1}_{2}'.format(output_directory, url, file_extension) for url in assemblies_ids]

    print('\nStarting download of {0} {1} files...'.format(files_number, file_extension))

    start: float = time.time()
    # We can use a with statement to ensure threads are cleaned up promptly
    failures: List[str] = []
    success: int = 0
    with concurrent.futures.ThreadPoolExecutor(max_workers=threads) as executor:
        # Start the load operations and mark each future with its URL
        for res in executor.map(df.download_file, ftp_urls, assemblies_ids, repeat(retry)):
            if 'Failed' in res:
                failures.append(res)
            else:
                success += 1
                print('\r', 'Downloaded {0}/{1}'.format(success, files_number), end='')

    print('\nFailed download for {0} files.'.format(len(failures)))

    end: float = time.time()
    delta: float = end - start
    minutes: int = int(delta / 60)
    seconds: float = delta % 60
    print('\nFinished downloading {0}/{1} fasta.gz files.'
          '\nElapsed Time: {2}m{3:.0f}s'.format(files_number - len(failures), files_number, minutes, seconds))
    
    with open(os.path.join(output_directory, "assemblies_ids_ncbi.tsv"), 'w+') as ids_to_tsv:
        ids_to_tsv.write("\n".join(map(str, assemblies_ids)))


def parse_arguments():

    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument('-t', '--input-table', type=str,
                        required=True, dest='input_table',
                        help='TSV file downloaded from the '
                             'NCBI Genome Assembly and Annotation '
                             'report.')

    parser.add_argument('-o', '--output-directory', type=str,
                        required=True, dest='output_directory',
                        help='Path to the directory to which '
                             'downloaded files will be stored.')

    parser.add_argument('-s', '--species', type=str,
                        required=True, dest='species',
                        help='Organism name. For example: '
                             'Streptococcus pneumoniae.')

    parser.add_argument('--fe', '--file-extension', type=str,
                        required=False,
                        choices=['genomic.fna.gz', 'assembly_report.txt',
                                 'assembly_status.txt', 'cds_from_genomic.fna.gz',
                                 'feature_count.txt.gz', 'feature_table.txt.gz',
                                 'genomic.gbff.gz', 'genomic.gff.gz',
                                 'genomic.gtf.gz', 'protein.faa.gz',
                                 'protein.gpff.gz', 'rna_from_genomic.fna.gz',
                                 'translated_cds.faa.gz'],
                        default='genomic.fna.gz',
                        dest='file_extension',
                        help='Choose file type to download through '
                             'extension.')

    parser.add_argument('--ftp', type=str, required=False,
                        choices=['refseq+genbank', 'refseq', 'genbank'],
                        default='refseq+genbank', dest='ftp',
                        help='The script can search for the files to '
                             'download in RefSeq or Genbank or both '
                             '(will only search in Genbank if download '
                             'from RefSeq fails).')

    parser.add_argument('-th', '--threads', type=int,
                        required=False, default=2,
                        dest='threads',
                        help='Number of threads for download.')

    parser.add_argument('-r', '--retry', type=int,
                        required=False, dest='retry',
                        default=7,
                        help='Maximum number of retries when a '
                             'download fails.')

    args = parser.parse_args()

    return args


if __name__ == '__main__':

    args = parse_arguments()
    main(**vars(args))
