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


def print_classifications_results(merged_all_classes: tp.MergedAllClasses, drop_possible_loci: List[int], 
                                  to_blast_paths: Dict[str, str], clusters: Dict[str, Any], count_classes_final: Dict[str, int]) -> None:
    """
    Prints the classification results based on the provided parameters.

    Parameters
    ----------
    merged_all_classes : tp.MergedAllClasses
        Dictionary containing CDS to keep, classified by their class type.
    drop_possible_loci : List[int]
        List of possible loci dropped.
    to_blast_paths : Dict[str, str]
        Path to BLAST.
    clusters : Dict[str, Any]
        The dictionary containing the clusters.
    moved_recs : Dict[str, List[Set[str]]]
        The dictionary of which loci in each class got the recommendation they were expected to or were moved, for exemple 'Drop' instead of 'Choice'.

    Returns
    -------
    None
        Prints the classification results to stdout.

    Notes
    -----
    - The function first processes the `clusters_to_keep` dictionary to count the number of CDS
      representatives for each class.
    - It then prints the results, including the number of groups classified under each class and
      any recommendations for verification.
    - If there are any retained groups not matched by BLASTn, it handles them separately.
    """
    def print_results(class_: str, count: int, printout: Dict[str, Any]) -> None:
        """
        Prints the classification results based on the class type.

        Parameters
        ----------
        class_ : str
            The class type.
        count : int
            The count of groups.
        printout : Dict[str, Any]
            The dictionary containing printout information.
        moved_recs : Dict[str, List[Set[str]]]
            The dictionary of which loci in each class got the recommendation they were expected to or were moved, for exemple 'Drop' instead of 'Choice'.

        Returns
        -------
        None
            Prints the classification results to stdout.

        Notes
        -----
        - The function prints different messages based on the class type.
        - It provides recommendations for verification for certain classes.
        """
        if count > 0:
            if class_ in ['2b', '4b']:
                prf.print_message(f"\t\t{count} loci are classified as {class_} and were retained"
                            " but it is recommended to verify them as they may be contained or contain partially inside"
                            " their BLAST match.", None)
            elif class_ == '1a':
                prf.print_message(f"\t\t{count} loci are classified as {class_}"
                            f" and are contained in {len(printout['1a'])} joined groups that were retained.", None)
            elif class_ in ['4c', '5']:
                prf.print_message(f"\t\t{count} loci are classified as {class_}. These will be added to recommendations with 'Add'.", None)
            else:
                prf.print_message(f"\t\t{count} loci are classified as {class_}.", None)

    # If 'Retained_not_matched_by_blastn' exists in clusters_to_keep, remove it and store it separately
    retained_not_matched_by_blastn: Optional[Any] = merged_all_classes.pop('Retained_not_matched_by_blastn', None)

    # Display info about the results obtained from processing the classes.
    # Get the total number of CDS reps considered for classification.
    count_cases: Dict[str, int] = {}
    # Check if loci is not empty
    total_loci: int = sum(count_cases.values()) + sum(count_classes_final.values())
    prf.print_message(f"Out of {len(to_blast_paths)}:", None)
    prf.print_message(f"\t{total_loci} representatives had matches with BLASTn against the schema DNA sequences.", None)
    for class_, count in count_classes_final.items():
        print_results(class_, count, merged_all_classes)
    for class_, count in count_cases.items():
        print_results(class_, count, merged_all_classes)
    prf.print_message(f"\tOut of those {len(to_blast_paths.values()) - sum(count_cases.values()) - sum(count_classes_final.values())} didn't have any matches", None)

    if retained_not_matched_by_blastn:
        merged_all_classes['Retained_not_matched_by_blastn'] = retained_not_matched_by_blastn


def prepare_loci(schema_folder: str,
                 constants: List[Any],
                 results_output: str) -> Tuple[
                     Dict[str, str], 
                     str, 
                     Dict[str, str],
                     Dict[str, str],  
                     Dict[str, List[str]],  
                     Dict[str, List[str]], 
                     Dict[str, List[str]],
                     Dict[str, str],
                     Dict[str, int],
                     Dict[str, int]]:
    """
    Process new loci by translating sequences, counting frequencies, and preparing files for BLAST.

    Parameters
    ----------
    schema_folder : str
        Path to the folder containing schema FASTA files.
    constants : list
        A list of constants used for processing.
    results_output : str
        Path to the directory where results will be saved.

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
