import os
import csv
import socket
import concurrent.futures
from itertools import repeat
from typing import List, Optional

try:
    from utils.download_functions import download_file
    from utils.file_functions import create_directory
    from utils.constants import SOCKET_TIMEOUT, PROTEOME_TEMPLATE_URL
    from utils import print_functions as pf
except ModuleNotFoundError:
    from SchemaRefinery.utils.download_functions import download_file
    from SchemaRefinery.utils.file_functions import create_directory
    from SchemaRefinery.utils.constants import SOCKET_TIMEOUT, PROTEOME_TEMPLATE_URL
    from SchemaRefinery.utils import print_functions as pf

# Set socket timeout for urllib calls
socket.setdefaulttimeout(SOCKET_TIMEOUT)

def proteome_fetcher(proteome_table: str, output_directory: str, threads: int, retry: int) -> Optional[str]:
    """
    Fetch proteomes from UniProt based on a table of proteome identifiers.

    Parameters
    ----------
    proteome_table : str
        Path to the input table containing proteome identifiers.
    output_directory : str
        Path to the output directory where proteomes will be saved.
    threads : int
        Number of threads to use for downloading.
    retry : int
        Number of retry attempts for failed downloads.

    Returns
    -------
    Optional[str]
        Path to the directory containing downloaded proteomes, or None if no proteomes were downloaded.
    """
    # Create directory for proteomes
    proteomes_directory: str = os.path.join(output_directory, 'proteomes')
    create_directory(proteomes_directory)

    # Open table downloaded from UniProt
    with open(proteome_table, 'r') as infile:
        reader = csv.reader(infile, delimiter='\t')
        # Exclude header
        next(reader, None)
        # Get proteome identifiers
        proteome_ids: List[str] = [row[0] for row in reader]

    # Build the UniProt URLs
    urls: List[str] = [PROTEOME_TEMPLATE_URL.format(pid) for pid in proteome_ids]
    num_proteomes: int = len(urls)
    pf.print_message(f'Proteomes to download: {num_proteomes}', 'info')

    if num_proteomes == 0:
        pf.print_message('Input proteome table contained no proteome ids.', 'warning')
        return None

    filenames: List[str] = ['{0}.fasta.gz'.format(pid) for pid in proteome_ids]
    filepaths: List[str] = [os.path.join(proteomes_directory, filename) for filename in filenames]

    pf.print_message(f"Starting download of {num_proteomes} proteomes...", 'info')

    # We can use a with statement to ensure threads are cleaned up promptly
    success: int = 0
    failures: List[str] = []
    with concurrent.futures.ThreadPoolExecutor(max_workers=threads) as executor:
        # Start the load operations and mark each future with its URL
        for res in executor.map(download_file, urls, filepaths, repeat(retry)):
            if 'Failed' in res:
                failures.append(res[0])
            elif os.path.getsize(res[1]) == 0:
                pf.print_message(f"Error: The downloaded archive '{res[1]}' is empty.", 'error')
                os.remove(res[1])
            else:
                success += 1
                pf.print_message(f'Downloaded {success}/{num_proteomes}', "info", end='\r', flush=True)

    pf.print_message(f"Failed download for {len(failures)} proteomes.", 'warning')

    return proteomes_directory