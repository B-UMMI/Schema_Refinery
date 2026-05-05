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

def write_dropped_cds_to_file(dropped_cds: Dict[str, str], results_output: str) -> None:
    """
    Write dropped CDS to file.

    Parameters
    ----------
    dropped_cds : Dict[str, str]
        The dictionary containing the dropped CDSs and their reasons.
    results_output : str
        The path where the output will be written.

    Returns
    -------
    None
        Writes to file.
    """
    dropped_cds_output: str = os.path.join(results_output, 'dropped_cds.tsv')
    with open(dropped_cds_output, 'w') as dropped_cds_file:
        dropped_cds_file.write('CDS_ID\tReason_for_dropping\n')
        for cds_id, reason in dropped_cds.items():
            dropped_cds_file.write(f"{cds_id}\t{reason}\n")

def update_ids_and_save_changes(merged_all_classes: tp.MergedAllClasses, 
                                clusters: Dict[str, List[str]], 
                                cds_original_ids: Dict[str, List[str]], 
                                dropped_cds: Dict[str, str],
                                all_nucleotide_sequences: Dict[str, str], 
                                results_output: str) -> None:
    """
    Update the IDs based on clustering and joining operations and save the changes.

    This function iterates through each class and its corresponding group of CDS (Coding DNA Sequences) to keep,
    updates the IDs based on the provided clusters and the original to new ID mappings, and saves the final ID changes
    to a TSV (Tab-Separated Values) file in the specified output directory.

    Parameters
    ----------
    merged_all_classes : tp.MergedAllClasses
        A dictionary where each key is a class and each value is a list of CDS to keep.
    clusters : Dict[str, List[str]]
        A dictionary mapping representative IDs to their cluster members.
    cds_original_ids : Dict[str, List[str]]
        A dictionary mapping original IDs to their new IDs after processing.
    dropped_cds : Dict[str, str]
        A dictionary mapping all of the dropped CDSs to the cause of drop.
    all_nucleotide_sequences : Dict[str, str]
        A dictionary that contains DNA sequences for each CDS.
    results_output : str
        The directory path where the ID changes file will be saved.

    Returns
    -------
    None
        The function writes the output files to the specified directory.

    Notes
    -----
    The function iterates through the `merged_all_classes` dictionary, updating IDs for each CDS based on their membership
    in the provided `clusters`. It generates a new ID for each CDS, updates `cds_original_ids` with these new IDs,
    and writes the original and new IDs to a TSV file named 'cds_id_changes.tsv' in the `results_output` directory.

    The ID updating process involves generating a new ID by appending an index to the main representative ID for each
    CDS in a cluster. This index is incremented for each CDS in the cluster.

    Examples
    --------
    Assuming the existence of appropriate dictionaries for `merged_all_classes`, `clusters`, `cds_original_ids`, and a valid
    path for `results_output`, the function can be called as follows:

    >>> update_ids_and_save_changes(merged_all_classes, clusters, cds_original_ids, dropped_cds, all_nucleotide_sequences, '/path/to/output')
    
    This would process the IDs as described and save the changes to '/path/to/output/cds_id_changes.tsv'.
    """

    # Iterate through each class and its CDS group
    for class_, cds_group in merged_all_classes.items():
        for cds in cds_group:
            main_rep: str = cds  # The main representative ID for the CDS group
            
            # If the class is not '1a', treat the CDS as a single-element list
            if class_ != '1a':
                cds = [cds]
            else:
                # For class '1a', get the CDS group from merged_all_classes
                cds = merged_all_classes[class_][cds]
            
            index: int = 1  # Initialize an index for creating new IDs
            
            # Iterate through each representative ID in the CDS group
            for rep_id in list(cds):
                # Get all CDS IDs in the cluster for the representative ID
                cds_ids: List[str] = clusters[rep_id]
                # Delete clusters with old IDs
                del clusters[rep_id]
                # Create new rep ID
                # Iterate through each CDS ID in the cluster
                for cds_id in list(cds_ids):
                    # Skip cases
                    if cds_id in dropped_cds:
                        continue
                    if not clusters.get(rep_id):
                        clusters[rep_id] = []
                    # Create a new ID using the main representative ID and the index
                    new_id: str = f"{main_rep}_{index}"
                    # Update the original ID with the new ID in cds_original_ids
                    cds_id_first: str = itf.identify_string_in_dict_get_key(cds_id, cds_original_ids)
                    cds_id_second: str = itf.identify_string_in_dict_get_value(cds_id, cds_original_ids)[-1]
                    # Replace in FASTA dict
                    all_nucleotide_sequences[new_id] = all_nucleotide_sequences.pop(cds_id_second)
                    # Add new cluster ID
                    clusters[rep_id].append(new_id)
                    cds_original_ids[cds_id_first].append(new_id)
                    index += 1  # Increment the index for the next ID
            
    # Prepare to write the ID changes to a file
    tab: str = "\t"
    id_changes_file: str = os.path.join(results_output, 'cds_id_changes.tsv')
    
    # Open the file and write the header and ID changes
    with open(id_changes_file, 'w') as id_changes:
        id_changes.write('Original_ID\tID_after_clustering\tID_after_joining\n')
        for original_ids, changed_ids in cds_original_ids.items():
            # Write each original ID and its changed IDs to the file
            id_changes.write(f"{original_ids}\t{tab.join(changed_ids)}\n")


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


def set_minimum_genomes_threshold(allelecall_directory: str, constants: List[Union[int, float]]) -> None:
    """
    Sets the minimum genomes threshold based on the dataset size.

    Parameters
    ----------
    allelecall_directory : str
        Path to the directory containing the allele calling results.
    constants : List[Union[int, float]]
        List of constants where the threshold will be set.

    Returns
    -------
    None
        The function updates the constants list in place.
    """
    try:
        genome_list: List[str] = pd.read_csv(ff.join_paths(allelecall_directory, ["results_statistics.tsv"]), sep='\t', usecols=['FILE'])['FILE'].tolist()
        number_of_genomes: int = len(genome_list)

        if number_of_genomes <= 20:
            constants[2] = 5
        else:
            constants[2] = round(number_of_genomes * 0.01)
    except Exception as e:
        pf.print_message(f"Setting minimum genomes threshold: {e}", 'error')
        constants[2] = 5  # Default value in case of error


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


def count_cds_frequency(unclassified_cds: Dict[str, str], decoded_seqids: Dict[str, List[float]]) -> Tuple[Dict[str, int], Dict[str, List[str]]]:
    """
    Counts the frequency of CDS in the genomes and identifies their presence.

    Parameters
    ----------
    unclassified_cds : Dict[str, str]
        Dictionary with CDS IDs as keys and sequences as values.
    decoded_seqids : Dict[str, List[str]]
        Dictionary with hashed sequences as keys and genome IDs as values.

    Returns
    -------
    Dict[str, int]
        A tuple containing the frequency of CDS and their presence in genomes.
    """
    cds_frequency: Dict[str, int] = {}
    # The unclassified CDSs and the decoded sequence IDs come from the same run
    # So each unclassified CDS should have a corresponding hashed sequence in the decoded sequence IDs
    for seqid, sequence in unclassified_cds.items():
        hashed_seq: str = sf.seq_to_hash(str(sequence))
        # Count only the unique genome IDs for the frequency [1:]
        cds_frequency[seqid] = len(set(decoded_seqids[hashed_seq][1:]))

    return cds_frequency


def remove_dropped_cds(all_translation_dict: Dict[str, str], 
                       dropped_cds: Dict[str, str], 
                       protein_hashes: Dict[str, List[str]]) -> Dict[str, str]:
    """
    Removes dropped CDS from the translation dictionary and updates protein hashes.

    Parameters
    ----------
    all_translation_dict : Dict[str, str]
        Dictionary with CDS IDs as keys and sequences as values.
    dropped_cds : Dict[str, str]
        Set of CDS IDs that are dropped.
    protein_hashes : Dict[str, List[str]]
        Dictionary with protein hashes as keys and CDS IDs as values.

    Returns
    -------
    Dict[str, str]
        Updated translation dictionary.
    """
    for key in list(all_translation_dict.keys()):
        if key in dropped_cds:
            protein_hash: str = itf.identify_string_in_dict_get_key(key, protein_hashes)
            same_protein_id: List[str] = protein_hashes[protein_hash]
            if key == same_protein_id[0]:
                protein_hashes[protein_hash].remove(key)
                if not protein_hashes[protein_hash]:
                    del protein_hashes[protein_hash]
                    continue
                new_id: str = protein_hashes[protein_hash][0]
                all_translation_dict[new_id] = all_translation_dict.pop(key)
            else:
                protein_hashes[protein_hash].remove(key)

    return {k: v for k, v in all_translation_dict.items() if k not in dropped_cds}


def get_representative_translation_dict(all_translation_dict: Dict[str, str], 
                                        all_alleles: Dict[str, List[str]]) -> Dict[str, str]:
    """
    Filters and sorts the representative translation dictionary.

    Parameters
    ----------
    all_translation_dict : Dict[str, str]
        Dictionary with CDS IDs as keys and sequences as values.
    all_alleles : Dict[str, List[str]]
        Dictionary with cluster representatives as keys and cluster members as values.

    Returns
    -------
    Dict[str, str]
        Sorted representative translation dictionary.
    """
    # Filter the representatives protein sequence
    reps_translation_dict: Dict[str, str] = {
        rep_id: rep_seq for rep_id, rep_seq in all_translation_dict.items()
        if rep_id.split('_')[0] in all_alleles
    }
    
    # Sort the representative translation dict from largest to smallest
    return {k: v for k, v in sorted(reps_translation_dict.items(), key=lambda x: len(x[1]), reverse=True)}


def calculate_kmers_similarity(sequences: Dict[str, str], 
                               representative_kmers: Dict[str, List[str]], 
                               sequence_lengths: Dict[str, int]) -> Dict[str, Dict[str, List[int]]]:
    """
    Calculates the k-mers similarity and coverage between representatives.

    Parameters
    ----------
    sequences : Dict[str, str]
        Dictionary with sequence IDs as keys and sequences as values.
    representative_kmers : Dict[str, List[str]]
        Dictionary with representative k-mers.
    sequence_lengths : Dict[str, int]
        Dictionary with sequence lengths.

    Returns
    -------
    Dict[str, Dict[str, List[int]]]
        Dictionary with k-mers similarity and coverage between representatives.
    """
    similar_sequences: Dict[str, Dict[str, List[int]]] = {}
    for seqid, sequence in sequences.items():
        minimizers: Set[str] = set(kf.determine_minimizers(sequence, 5, 5, 1, True, True))
        similar_sequences[seqid] = cf.select_representatives(minimizers, representative_kmers, 0, 0, sequence_lengths, seqid, 5, False)

    return similar_sequences


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


def prepare_files_to_blast(output_directory: str,
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


def update_frequencies_in_genomes(clusters_to_keep_1a: Dict[str, List[str]], 
                                  frequency_in_genomes: Dict[str, int]) -> Dict[str, int]:
    """
    Updates the frequencies in genomes for joined groups and updates the changed clusters frequency from joined CDSs.

    Parameters
    ----------
    clusters_to_keep : Dict[str, List[str]]
        Dictionary containing clusters to keep with their members.
    frequency_in_genomes : Dict[str, int]
        Dictionary with the frequency of CDS in the genomes.

    Returns
    -------
    Dict[str, int]
        Updated frequency of CDS in the genomes.
    """
    updated_frequency_in_genomes: Dict[str, int] = {}
    new_cluster_freq: Dict[str, int] = {}

    # Calculate new frequencies for joined groups
    for cluster_id, cluster_members in clusters_to_keep_1a.items():
        new_cluster_freq[cluster_id] = sum(frequency_in_genomes[member] for member in cluster_members)
        for member in cluster_members:
            updated_frequency_in_genomes[member] = new_cluster_freq[cluster_id]

    # Add all the other frequencies
    updated_frequency_in_genomes.update(frequency_in_genomes)
    updated_frequency_in_genomes.update(new_cluster_freq)

    return updated_frequency_in_genomes
