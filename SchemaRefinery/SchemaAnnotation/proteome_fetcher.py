#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Purpose
-------

This script accepts a TSV file downloaded from UniProt with
information about a set of proteomes and downloads the
proteomes listed in the file.

Code documentation
------------------
"""


import os
import csv
import socket
import concurrent.futures
from itertools import repeat

try:
    from utils.download_functions import download_file
    from utils.file_functions import create_directory
    from utils.constants import SOCKET_TIMEOUT, PROTEOME_TEMPLATE_URL
except ModuleNotFoundError:
    from SchemaRefinery.utils.download_functions import download_file
    from SchemaRefinery.utils.file_functions import create_directory
    from SchemaRefinery.utils.constants import SOCKET_TIMEOUT, PROTEOME_TEMPLATE_URL


# Set socket timeout for urllib calls
socket.setdefaulttimeout(SOCKET_TIMEOUT)


def proteome_fetcher(proteome_table: str, output_directory: str, threads: int,
                     retry: int):

    proteomes_directory = os.path.join(output_directory, 'proteomes')
    create_directory(proteomes_directory)

    # Open table downloaded from UniProt
    with open(proteome_table, 'r') as infile:
        reader = csv.reader(infile, delimiter='\t')
        # Exclude header
        next(reader, None)
        # Get proteome identifiers
        proteome_ids = [row[0] for row in reader]

    # Build the Uniprot URLs
    urls = [PROTEOME_TEMPLATE_URL.format(pid) for pid in proteome_ids]
    num_proteomes = len(urls)
    print(f'Proteomes to download: {num_proteomes}')

    if num_proteomes == 0:
        print('Input proteome table contained no proteome ids.')
        return

    filenames = ['{0}.fasta.gz'.format(pid) for pid in proteome_ids]
    filepaths = [os.path.join(proteomes_directory, filename)
                 for filename in filenames]

    print('\nStarting download of {0} proteomes...'.format(num_proteomes))

    # We can use a with statement to ensure threads are cleaned up promptly
    success = 0
    failures = []
    with concurrent.futures.ThreadPoolExecutor(max_workers=threads) as executor:
        # Start the load operations and mark each future with its URL
        for res in executor.map(download_file, urls, filepaths, repeat(retry)):
            if 'Failed' in res:
                failures.append(res)
            else:
                success += 1
                print('\r', 'Downloaded {0}/{1}'.format(success, num_proteomes), end='')

    print(f'\nFailed download for {len(failures)} proteomes.')

    return proteomes_directory
