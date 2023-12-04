import os
import sys
import concurrent.futures
from itertools import repeat

try:
    from utils.sequence_functions import read_fasta_file, seq_to_hash
except ModuleNotFoundError:
    from SchemaRefinery.utils.sequence_functions import read_fasta_file, seq_to_hash

def hash_sequences(file_paths):
    """
    Hashes sequences in fasta files based on input file_paths.

    Parameters
    ----------
    file_paths : list
        Contains file path to the fasta files.

    Returns
    -------
    return : list
        Returns a list containing all of the sequences hashes present in the input files.
    """

    hash_list = []
    for fasta in file_paths.values():

        for seq in read_fasta_file(fasta).values():
            hash_list += seq_to_hash(seq)

    return hash_list

def not_included_cds(hash_list, file_path_cds):
    """
    Compares the hashes list with the hashes obtained from cds.

    Parameters
    ----------
    hash_list : list
        List containing all of the sequences hashes present in the input files.

    Returns
    -------
    return : dict
        Returns dict with key as fasta header and value as fasta sequence.
    """

    not_included_cds = {}
    for record in read_fasta_file(file_path_cds):
        if seq_to_hash(str(record.seq)) not in hash_list:
            not_included_cds[record.id] = str(record.seq)

    return not_included_cds

def main(schema, output_directory, allelecall_directory):

    schema_file_paths = {f.replace(".fasta", ""): os.path.join(schema, f) for f in os.listdir(schema) if not os.path.isdir(f) and f.endswith(".fasta")}
    file_path_cds = os.path.join(allelecall_directory, "temp", "3_cds_preprocess", "cds_deduplication", "distinct_cds_merged.fasta")

    if not os.path.exists(file_path_cds):
        sys.exit(f"Error: {file_path_cds} must exist, make sure that allele call was run using --no-cleanup flag.")

    not_included_cds = not_included_cds(hash_sequences(schema_file_paths), file_path_cds)
