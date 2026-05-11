import os
import time
import shutil
from typing import Dict, Any, List, Tuple, Union, Optional, Set

import pandas as pd

try:
    from utils import (file_functions as ff,
                                        clustering_functions as cf,
                                        blast_functions as bf,
                                        alignments_functions as af,
                                        iterable_functions as itf,
                                        linux_functions as lf,
                                        classify_cds_functions as ccf,
                                        graphical_functions as gf,
                                        pandas_functions as pf,
                                        sequence_functions as sf,
                                        Types as tp,
                                        print_functions as prf)
except ModuleNotFoundError:
    from SchemaRefinery.utils import (file_functions as ff,
                                        clustering_functions as cf,
                                        blast_functions as bf,
                                        alignments_functions as af,
                                        iterable_functions as itf,
                                        linux_functions as lf,
                                        classify_cds_functions as ccf,
                                        graphical_functions as gf,
                                        pandas_functions as pf,
                                        sequence_functions as sf,
                                        Types as tp,
                                        print_functions as prf)


def prepare_cluster_blast_infiles(output_directory: str,
                           clusters: Dict[str, List[str]], 
                           sequences: Dict[str, str],
                           protein_sequences: Dict[str, str],
                           blastdb_aliastool_exec) -> Dict[str, str]:
    """
    Writes representative and master FASTA files for BLAST.

    Parameters
    ----------
    output_directory : str
        Path to the directory where the BLAST input files will be stored.
    clusters : Dict[str, List[str]]
        Dictionary with cluster representatives as keys and cluster members as values.
    sequences : Dict[str, str]
        Dictionary with sequence IDs as keys and sequences as values.
    protein_sequences : Dict[str, str]
        Dictionary with protein sequence IDs as keys and protein sequences as values.

    Returns
    -------
    file_paths : Dict[str, List[str]]
        Dictionary containing repid as keys and a list of paths to the corresponding DNA FASTA, protein FASTA, and seqid TXT files as values.
    """
    file_paths: Dict[str, List[str]] = {}
    for repid, members in clusters.items():
        dna_fasta: str = os.path.join(output_directory, f"{repid}.fasta")
        protein_fasta: str = os.path.join(output_directory, f"{repid}_protein.fasta")
        seqid_file: str = os.path.join(output_directory, f"{repid}.txt")
        seqid_file_alias: str = os.path.join(output_directory, f"{repid}_alias")
        file_paths[repid] = [dna_fasta, protein_fasta, seqid_file, seqid_file_alias]

        # Create FASTA file with DNA sequence of the representative
        sequence: Dict[str, str] = {repid: sequences[repid]}
        ccf.write_fasta_file(dna_fasta, sequence)

        # Create FASTA file with protein sequence of the representative
        protein: Dict[str, str] = {repid: protein_sequences[repid]}
        ccf.write_fasta_file(protein_fasta, protein)

        # Create TXT file with the representative seqid
        with open(seqid_file, 'w') as outfile:
            outfile.write(f"{repid}\n")

        # Run blastdb_aliastool to create the alias files for the seqid files
        stdout, stderr = bf.run_blastdb_aliastool_multiprocessing(blastdb_aliastool_exec, seqid_file, seqid_file_alias)

    # Concatenate all the representative FASTA files into one file for BLASTn
    concatenated_dna_fasta: str = os.path.join(output_directory, "concatenated_dna.fasta")
    ff.concatenate_files([file[0] for file in file_paths.values()], concatenated_dna_fasta, header=None)

    # Concatenate all the representative FASTA files into one file for BLASTp
    concatenated_protein_fasta: str = os.path.join(output_directory, "concatenated_protein.fasta")
    ff.concatenate_files([file[1] for file in file_paths.values()], concatenated_protein_fasta, header=None)

    return file_paths, (concatenated_dna_fasta, concatenated_protein_fasta)


def prepare_loci_blast_infiles(schema_folder: str, output_directory: str, minimum_length: int, translation_table: int, blastdb_aliastool_exec: str):
    """
    Process new loci by translating sequences, counting frequencies, and preparing files for BLAST.

    Parameters
    ----------
    schema_folder : str
        Path to the folder containing schema FASTA files.
    output_directory : str
        Path to the directory where results will be saved.
    constants : list
        A list of constants used for processing.

    Returns
    -------
    tuple
        A tuple containing:
        - all_nucleotide_sequences (Dict[str, str]): Dictionary of nucleotide sequences.
        - master_file_path (str): Path to the master FASTA file.
        - trans_paths (Dict[str, str]): Dictionary with the path of the translation file for each loci.
        - to_blast_paths (Dict[str, str]): Dictionary of paths to sequences to be used for BLAST.
        - all_alleles (Dict[str, List[str]]): Dictionary of all alleles with loci IDs as keys.
        - group_reps_ids (Dict[str, List[str]]): Dictionary of group representative IDs.
        - group_alleles_ids (Dict[str, List[str]]): Dictionary of group allele IDs.
        - to_run_against (Dict[str, str]): Dictionary of paths to sequences to be used for BLAST database.
        - new_max_hits (Dict[str, int]): Dictionary with the number of max hits per loci for the BLAST.
        - seqid_file_dict (Dict[str, int]): Dictionary with the paths to the files with the seqid to ignore per loci for the BLAST.
    """
    # Map loci IDs to the paths to the FASTA files in the schema
    loci_files = {ff.file_basename(file, False): os.path.join(schema_folder, file) for file in os.listdir(schema_folder) if file.endswith('.fasta')}
    # Do the same for the FASTA files containing the representative alleles
    schema_short_dir = os.path.join(schema_folder, 'short')
    representative_files = {file.rsplit('_', 1)[0]: os.path.join(schema_short_dir, file) for file in os.listdir(schema_short_dir) if file.endswith('.fasta')}
    file_paths: Dict[str, List[str]] = {}
    for locus, file in representative_files.items():
        dna_fasta = file
        protein_fasta: str = os.path.join(output_directory, f"{locus}_protein.fasta")
        seqid_file: str = os.path.join(output_directory, f"{locus}.txt")
        seqid_file_alias: str = os.path.join(output_directory, f"{locus}_alias")
        file_paths[locus] = [dna_fasta, protein_fasta, seqid_file, seqid_file_alias]

        # Create FASTA file with protein sequence of the representative
        representative_alleles = sf.read_fasta_file_dict(dna_fasta, False)
        translated_alleles, _ = sf.translate_seq_deduplicate(representative_alleles, output_directory, minimum_length, translation_table, False)
        ccf.write_fasta_file(protein_fasta, translated_alleles)

        # Create TXT file with the representative seqid
        representative_seqids = '\n'.join(list(representative_alleles.keys()))
        with open(seqid_file, 'w') as outfile:
            outfile.write(f"{representative_seqids}\n")

        # Run blastdb_aliastool to create the alias files for the seqid files
        stdout, stderr = bf.run_blastdb_aliastool_multiprocessing(blastdb_aliastool_exec, seqid_file, seqid_file_alias)

    # Concatenate all the representative FASTA files into one file for BLASTn
    concatenated_dna_fasta: str = os.path.join(output_directory, "concatenated_dna.fasta")
    ff.concatenate_files([file[0] for file in file_paths.values()], concatenated_dna_fasta, header=None)

    # Concatenate all the representative FASTA files into one file for BLASTp
    concatenated_protein_fasta: str = os.path.join(output_directory, "concatenated_protein.fasta")
    ff.concatenate_files([file[1] for file in file_paths.values()], concatenated_protein_fasta, header=None)

    return file_paths, (concatenated_dna_fasta, concatenated_protein_fasta)
