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
        write_fasta_file(dna_fasta, sequence)

        # Create FASTA file with protein sequence of the representative
        protein: Dict[str, str] = {repid: protein_sequences[repid]}
        write_fasta_file(protein_fasta, protein)

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


def prepare_loci_blast_infiles(schema_folder: str, output_directory: str, constants: List[Any]):
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
    loci_files = {file.rsplit('_', 1)[0]: os.path.join(schema_folder, file) for file in os.listdir(schema_folder) if file.endswith('.fasta')}
    # Do the same for the FASTA files containing the representative alleles
    schema_short_dir = os.path.join(schema_folder, 'short')
    representative_files = {file.rsplit('_', 1)[0]: os.path.join(schema_short_dir, file) for file in os.listdir(schema_short_dir) if file.endswith('.fasta')}

    # Path to the master FASTA file
    master_file_path = os.path.join(results_output, 'master_protein.fasta')
    seqid_output = os.path.join(results_output, 'seqid_files')
    ff.create_directory(seqid_output)

    # Create a directory for possible new loci translations
    possible_new_loci_translation_folder = os.path.join(results_output, 'schema_translation_folder')
    ff.create_directory(possible_new_loci_translation_folder)

    # Initialize dictionaries for alleles, translations, and frequencies
    all_alleles: Dict[str, List[str]] = {}
    all_nucleotide_sequences: Dict[str, str] = {} 
    trans_paths: Dict[str, str] = {}
    group_reps_ids: Dict[str, List[str]] = {} 
    group_alleles_ids: Dict[str, List[str]] = {} 
    blast_alleles: Dict[str, List[str]] = {}
    seqid_file_dict: Dict[str, str] = {}

    # Process alleles to run, DNA sequences
    for locus, file in representative_files.items():
        sequences = sf.fetch_fasta_dict(file, False)
        for seqid, sequence in sequences.items():
            group_reps_ids.setdefault(locus, []).append(seqid)
            all_nucleotide_sequences.setdefault(locus, str(sequence))

    # Write master file to run against, DNA sequences
    prf.print_message('Write master file for Blastp.', 'info')

    for loci, loci_path in to_run_against.items():
        negative_seqid_file = os.path.join(seqid_output, f'negative_seqid_{loci_id}.txt')
        blast_alleles.setdefault(loci_id, [])
        seqid_file_dict.setdefault(loci_id, negative_seqid_file)
        fasta_dict = sf.fetch_fasta_dict(loci_path, False)
        for allele_id, sequence in fasta_dict.items():
            blast_alleles[loci_id].append(allele_id)
            group_alleles_ids.setdefault(loci_id, []).append(allele_id)
            all_nucleotide_sequences.setdefault(allele_id, str(sequence))
            protseq = sf.translate_sequence(str(sequence), constants[6])
            # Write to master file
            write_type = 'a' if os.path.exists(master_file_path) else 'w'
            with open(master_file_path, write_type) as master_file:
                master_file.write(f">{allele_id}\n{str(protseq)}\n")
            write_type2 = 'a' if os.path.exists(negative_seqid_file) else 'w'
            with open(negative_seqid_file, write_type2) as seqid_file:
                seqid_file.write(f"{allele_id}\n")

    for loci, loci_path in to_blast_paths.items():
        loci_id = ff.file_basename(loci, False)
        all_alleles.setdefault(loci_id, [])
        fasta_dict = sf.fetch_fasta_dict(loci_path, False)
        for allele_id, sequence in fasta_dict.items():  
            all_alleles[loci_id].append(allele_id)
        # Translate sequences and update translation dictionary
        trans_path_file = os.path.join(possible_new_loci_translation_folder, f"{loci_id}.fasta")
        trans_paths[loci_id] = trans_path_file
        trans_dict, _, _ = sf.translate_seq_deduplicate(fasta_dict,
                                                        trans_path_file,
                                                        constants[5],
                                                        constants[6],
                                                        False)

    return all_nucleotide_sequences, master_file_path, trans_paths, to_blast_paths, all_alleles, group_reps_ids, group_alleles_ids, to_run_against, seqid_file_dict
