import os
from typing import Dict, List, Set, Union, Tuple

import pandas as pd

try:
    from utils import (file_functions as ff,
                       blast_functions as bf,
                       sequence_functions as sf,
                       clustering_functions as cf,
                       iterable_functions as itf,
                       kmers_functions as kf,
                       Types as tp,
                       print_functions as pf)
except ModuleNotFoundError:
    from SchemaRefinery.utils import (file_functions as ff,
                                      blast_functions as bf,
                                      sequence_functions as sf,
                                      clustering_functions as cf,
                                      iterable_functions as itf,
                                      kmers_functions as kf,
                                      Types as tp,
                                      print_functions as pf)


def write_fastas_to_files(clusters: Dict[str, List[str]], all_nucleotide_sequences: Dict[str, str], output_path: str) -> str:
    """
    Write the cluster sequences to FASTA files.

    Parameters
    ----------
    clusters : Dict[str, List[str]]
        Dictionary of clusters with their members.
    all_nucleotide_sequences : Dict[str, str]
        Dictionary of allele IDs and their DNA sequences.
    output_path : str
        The directory path where the output files will be saved.
    
    Returns
    -------
    temp_fastas_folder : str
        The path to the folder containing the FASTA files.
    """

    # Create a directory for the cluster FASTA files
    temp_fastas_folder: str = os.path.join(output_path, 'temp_fastas')
    ff.create_directory(temp_fastas_folder)

    # Write the cluster sequences to FASTA files
    for cluster, members in clusters.items():
        cluster_file: str = os.path.join(temp_fastas_folder, f'{cluster}.fasta')
        with open(cluster_file, 'w') as cluster_fasta:
            for member in members:
                cluster_fasta.write(f">{member}\n{all_nucleotide_sequences[member]}\n")
    
    return temp_fastas_folder


def set_minimum_genomes_threshold(allelecall_directory: str, user_frequency: int) -> None:
    """
    Determine the minimum CDS/allele frequency based on the dataset size.

    Parameters
    ----------
    allelecall_directory : str
        Path to the directory containing the allele calling results.

    Returns
    -------
    minimum_frequency : int
        The minimum frequency value for loci presence.
    """
    genome_list: List[str] = pd.read_csv(ff.join_paths(allelecall_directory, ["results_statistics.tsv"]), sep='\t', usecols=['FILE'])['FILE'].tolist()
    number_of_genomes: int = len(genome_list)

    if not user_frequency:
        # Define minimum frequency value (this will round to 0 if dataset has less than 26 genomes)
        minimum_frequency = round(number_of_genomes*0.02)
    else:
        if user_frequency > number_of_genomes:
            pf.print_message("Value passed for the minimum locus frequency in the genomes exceeds the size of the dataset. Please provide a smaller value.", "error")
            sys.exit(0)
        else:
            minimum_frequency = user_frequency

    return minimum_frequency


def filter_by_size(sequences: Dict[str, str], size_threshold: int) -> Tuple[Dict[str, str], Dict[str, str]]:
    """
    Filters CDS by their size and updates the dropped CDS dictionary.

    Parameters
    ----------
    sequences : Dict[str, str]
        Dictionary with CDS IDs as keys and sequences as values.
    size_threshold : int
        Minimum size threshold for CDS.

    Returns
    -------
    Tuple[Dict[str, int], Dict[str, str], Dict[str, str]]
        A tuple containing the CDS size dictionary, filtered CDS dictionary, and the dropped CDS dictionary.
    """
    total: int = len(sequences)
    filtered_sequences: Dict[str, str] = {seqid: sequence for seqid, sequence in sequences.items() if len(sequence) >= size_threshold}
    dropped: Set[str] = set(sequences.keys()) - set(filtered_sequences.keys())

    return filtered_sequences, dropped


def write_cds_to_fasta(all_nucleotide_sequences: Dict[str, str], output_path: str) -> None:
    """
    Writes CDS sequences to a FASTA file.

    Parameters
    ----------
    all_nucleotide_sequences : Dict[str, str]
        Dictionary with CDS IDs as keys and sequences as values.
    output_path : str
        Path to the output FASTA file.

    Returns
    -------
    None
        The function writes the sequences to the specified output file.
    """
    with open(output_path, 'w+') as cds_not_found:
        for id_, sequence in all_nucleotide_sequences.items():
            cds_not_found.write(f">{id_}\n{str(sequence)}\n")


def write_fasta_file(file_path: str, sequences: Dict[str, str]) -> None:
    """
    Writes sequences to a FASTA file.

    Parameters
    ----------
    file_path : str
        Path to the output FASTA file.
    sequences : Dict[str, str]
        Dictionary with sequence IDs as keys and sequences as values.

    Returns
    -------
    None
        The function writes the sequences to the specified output file.
    """
    write_type: str = 'a' if os.path.exists(file_path) else 'w'
    with open(file_path, write_type) as fasta_file:
        for seq_id, sequence in sequences.items():
            fasta_file.write(f">{seq_id}\n{sequence}\n")
