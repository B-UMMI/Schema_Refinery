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
import sys
import csv
import time
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


# set socket timeout for urllib calls
socket.setdefaulttimeout(SOCKET_TIMEOUT)

def proteome_fetcher(input_table:str, output_directory:str, threads:int, retry:int):

    output_directory = os.path.join(output_directory, "Downloaded_proteomes")
    create_directory(output_directory)

    # open table downloaded from Uniprot
    with open(input_table, 'r') as table:
        reader = csv.reader(table, delimiter=',')
        # exclude header
        next(reader, None)
        proteomes = [row[0].split("\t")[0] for row in reader]

    # Build the Uniprot URLs
    urls = [PROTEOME_TEMPLATE_URL.format(proteome_id)
            for proteome_id in proteomes]
    
    print("urls: ", urls)

    files_number = len(urls)
    if files_number == 0:
        sys.exit('No valid URLs.')

    local_filenames = ['{0}.fasta.gz'.format(proteome_id)
                       for proteome_id in proteomes]
    local_filepaths = [os.path.join(output_directory, filename)
                       for filename in local_filenames]

    print('\nStarting download of {0} proteomes...'.format(len(urls)))
    start = time.time()

    # We can use a with statement to ensure threads are cleaned up promptly
    failures = []
    success = 0
    with concurrent.futures.ThreadPoolExecutor(max_workers=threads) as executor:
        # Start the load operations and mark each future with its URL
        for res in executor.map(download_file, urls, local_filepaths, repeat(retry)):
            if 'Failed' in res:
                failures.append(res)
            else:
                success += 1
                print('\r', 'Downloaded {0}/{1}'.format(success, files_number), end='')

    print('\nFailed download for {0} files.'.format(len(failures)))

    end = time.time()
    delta = end - start
    minutes = int(delta/60)
    seconds = delta % 60
    print('\nFinished downloading {0}/{1} fasta.gz files.'
          '\nElapsed Time: {2}m{3:.0f}s'
          ''.format(files_number-len(failures),
                    files_number, minutes, seconds))
    
    return output_directory
