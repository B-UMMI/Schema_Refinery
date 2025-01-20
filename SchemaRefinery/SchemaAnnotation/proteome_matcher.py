import os
import pickle
from typing import Dict, List, Set, Tuple, Optional, Union

try:
    from utils import (sequence_functions as sf,
                       file_functions as ff,
                       clustering_functions as cf,
                       blast_functions as bf,
                       linux_functions as lf,
                       alignments_functions as af,
                       iterable_functions as itf,
                       pandas_functions as pf,
                       print_functions as prf)
except ModuleNotFoundError:
    from SchemaRefinery.utils import (sequence_functions as sf,
                                      file_functions as ff,
                                      clustering_functions as cf,
                                      blast_functions as bf,
                                      linux_functions as lf,
                                      alignments_functions as af,
                                      iterable_functions as itf,
                                      pandas_functions as pf,
                                      print_functions as prf)

def create_database_files(proteome_file: str, clustering_sim: float, clustering_cov: float, size_ratio: float,
                          blast_processing_folder: str) -> Union[tuple[str, Dict[str, List[str]], Dict[str, str]], Tuple[None, None, None]]:
    """
    Create database files from a proteome file by extracting protein sequences, clustering them, and creating a BLAST database.

    Parameters
    ----------
    proteome_file : str
        Path to the proteome file.
    clustering_sim : float
        Clustering similarity threshold.
    clustering_cov : float
        Clustering coverage threshold.
    size_ratio : float
        Size ratio for clustering.
    blast_processing_folder : str
        Path to the output directory.

    Returns
    -------
    Optional[str]
        Path to the created BLAST database files, or None if no proteins were found.
    """
    file_name: str = os.path.basename(proteome_file)
    
    # Initialize dictionaries and sets for storing protein sequences and annotations
    all_translation_dict: Dict[str, str] = {}
    processed_proteins: Set[str] = set()
    # Create dictionaries to hashed sequence and its representative ID
    hash_to_rep_id: Dict[str, str] = {}
    # Create dictionaries to store the representatives ID and what it representes
    same_protein_other_annotations: Dict[str, List[str]] = {}

    prf.print_message(f"Extracting protein from proteome file: {os.path.basename(proteome_file)}...", "info")
    
    # Fetch protein sequences from the proteome file
    fasta_dict: Dict[str, str] = sf.fetch_fasta_dict(proteome_file, False)
    total_proteins: int = len(fasta_dict)
    
    if total_proteins == 0:
        prf.print_message(f"No proteins found in {os.path.basename(proteome_file)}", "warning")
        return (None, None, None)

    # Process each protein sequence
    for id_, sequence in fasta_dict.items():
        protein_hash: str = sf.hash_sequence(sequence)
        if protein_hash in processed_proteins:
            # If the protein has already been processed, save the other annotations that it may have
            same_protein_other_annotations.setdefault(hash_to_rep_id[protein_hash], []).append(id_)
            continue
        else:
            processed_proteins.add(protein_hash)
            all_translation_dict[id_] = sequence
            hash_to_rep_id[protein_hash] = id_

    prf.print_message(f"Out of {total_proteins} protein sequences, {len(all_translation_dict)} are unique proteins", "info")
    
    prf.print_message(f"Clustering protein sequences...", "info")
    all_alleles: Dict[str, str] = {}
    reps_sequences: Dict[str, str] = {}
    reps_groups: Dict[str, List[str]] = {}
    prot_len_dict: Dict[str, int] = {}
    
    # Cluster protein sequences
    all_alleles, reps_sequences, reps_groups, prot_len_dict = cf.minimizer_clustering(
        all_translation_dict,
        5,
        5,
        True,
        1,
        all_alleles,
        reps_sequences,
        reps_groups,
        1,
        clustering_sim,
        clustering_cov,
        True,
        size_ratio
    )

    prf.print_message(f"Clustered {len(all_translation_dict)} into {len(reps_sequences)} clusters", "info")
    
    # Save clustered protein sequences to a file
    clustered_protein_master_file: str = os.path.join(blast_processing_folder, f"{file_name}")
    with open(clustered_protein_master_file, 'w') as outfile:
        for protein_id, values in reps_sequences.items():
            outfile.write(f">{protein_id}\n{values}\n")
    
    # Create BLAST database
    blastdb_path: str = os.path.join(blast_processing_folder, 'blastdb')
    ff.create_directory(blastdb_path)
    blast_db_files: str = os.path.join(blastdb_path, 'genbank_protein_db')
    makeblastdb_exec: str = lf.get_tool_path('makeblastdb')
    bf.make_blast_db(makeblastdb_exec, clustered_protein_master_file, blast_db_files, 'prot')
        
    return blast_db_files, same_protein_other_annotations, all_alleles

def run_blast_for_proteomes(max_id_length: int, proteome_file_ids: Dict[str, List[str]],
                            best_bsr_values_per_proteome_file: Dict[str, Dict[str, Tuple[str, float]]],
                            blast_processing_folder: str, translations_paths: Dict[str, str],
                            blast_db_files: Optional[str], self_score_dict: Dict[str, float], cpu: int,
                            bsr: float) -> dict[str, Tuple[str, float]]:
    """
    Run BLAST for proteomes to calculate self-scores and BSR values, and save the annotations.

    Parameters
    ----------
    max_id_length : int
        Maximum length of the ID values
    blast_processing_folder : str
        Path to the folder for BLAST processing.
    translations_paths : Dict[str, str]
        Dictionary of translation paths.
    blast_db_files : str
        Path to the BLAST database files.
    proteome_folder : str
        Path to the proteome folder.
    file_name_without_extension : str
        File name without extension.
    descriptions : Dict[str, str]
        Dictionary of descriptions.
    self_score_dict : Dict[str, float]
        Dictionary of self-scores for each representative ID.
    cpu : int
        Number of CPUs to use.
    bsr : float
        BSR threshold value.

    Returns
    -------
    None
    """
    # Get Path to the blastp executable
    get_blastp_exec: str = lf.get_tool_path('blastp')
    
    # Run BLASTp
    prf.print_message("Running BLASTp...", "info")
    blastp_results_folder: str = os.path.join(blast_processing_folder, 'blastp_results')
    ff.create_directory(blastp_results_folder)
    
    # Run BLASTp between all BLASTn matches (rep vs all its BLASTn matches).
    bsr_values: Dict[str, Dict[str, float]] = {}
    best_bsr_values: Dict[str, Tuple[str, float]] = {}
    total_blasts: int = len(translations_paths)
    # Run BLASTp in parallel
    blastp_results_files = bf.run_blastp_operations(cpu,
                                                    get_blastp_exec,
                                                    blast_db_files,
                                                    translations_paths,
                                                    blastp_results_folder,
                                                    total_blasts,
                                                    max_id_length)

    for blast_results_file in blastp_results_files:
        # Get the alignments
        filtered_alignments_dict: Dict[str, Dict[str, Dict[str, Dict[str, float]]]]
        filtered_alignments_dict, _, _, _ = af.get_alignments_dict_from_blast_results(blast_results_file, 0, True, False, True, True, False)

        # Since BLAST may find several local alignments, choose the largest one to calculate BSR.
        for query, subjects_dict in filtered_alignments_dict.items():
            # Get the loci name
            loci: str = query.split('_')[0]
            # Create the dict of the query
            bsr_values.setdefault(query, {})
            for subject_id, results in subjects_dict.items():
                # Highest score (First one)
                subject_score: float = next(iter(results.values()))['score']
                # Calculate BSR value
                bsr_value: float = bf.compute_bsr(subject_score, self_score_dict[query])
                # Check if the BSR value is higher than the threshold
                if bsr_value >= bsr:
                    # Round BSR values if they are superior to 1.0 to 1 decimal place
                    if bsr_value > 1.0:
                        bsr_value = round(bsr_value, 1)
                    # Save all of the different matches that this query had and their BSR values
                    bsr_values[query].update({subject_id: bsr_value})
                else:
                    continue
                # Check if the BSR value is the best for the locus
                current_best_bsr: Optional[Tuple[str, float]] = best_bsr_values.get(loci)
                # If there is a previous BSR value for the locus, check if the current BSR value is higher
                # We are interested in the best match only
                if not current_best_bsr:
                    best_bsr_values[loci] = (subject_id, bsr_value)
                elif bsr_value > current_best_bsr[1]:
                    best_bsr_values[loci] = (subject_id, bsr_value)
                    
                # Get best value for genbank file
                proteome_file: str = itf.identify_string_in_dict_get_key(subject_id, proteome_file_ids)
                current_best_in_proteome_file: Optional[Tuple[str, float]] = best_bsr_values_per_proteome_file[proteome_file].get(loci)
                if current_best_in_proteome_file and bsr_value > current_best_in_proteome_file[1]:
                    best_bsr_values_per_proteome_file[proteome_file][loci] = (subject_id, bsr_value)
                else:
                    best_bsr_values_per_proteome_file[proteome_file][loci] = (subject_id, bsr_value)

    return best_bsr_values

def proteome_matcher(proteome_files: List[str], proteome_file_ids: Dict[str, List[str]],
                    schema_directory: str, output_directory: str, cpu: int, bsr: float, 
                    translation_table: int, clustering_sim: float, 
                    clustering_cov: float, size_ratio: float, run_mode: str,
                    proteome_ids_to_add: List[str]) -> Tuple[str, str]:
    """
    Match proteomes by creating BLAST database files, translating sequences, and running BLAST.

    Parameters
    ----------
    proteome_files : List[str]
        List of paths to proteome files.
    schema_directory : str
        Path to the schema directory.
    output_directory : str
        Path to the output directory.
    cpu : int
        Number of CPUs to use.
    bsr : float
        BSR threshold value.
    translation_table : int
        Translation table number.
    clustering_sim : float
        Clustering similarity threshold.
    clustering_cov : float
        Clustering coverage threshold.
    size_ratio : float
        Size ratio for clustering.
    run_mode : str
        Mode to run ('alleles' or 'reps').

    Returns
    -------
    None
    """
    proteome_matcher_output = os.path.join(output_directory, 'proteome_matcher_output')
    # Create BLAST database files for each proteome
    proteomes_data_paths: Dict[str, List[Optional[str]]] = {}
    proteomes_data: List[Tuple[Optional[Dict[str, List[str]]], Optional[Dict[str, str]]]] = []
    for proteome_file in proteome_files[:2]:
        # Get proteome file name
        proteome_file_base: str = os.path.basename(proteome_file)
        # Create folder for proteome processing
        proteome_folder: str = os.path.join(proteome_matcher_output, f"{proteome_file_base.split('.')[0]}_processing")
        # Create directory for proteome BLAST processing
        blast_processing_folder: str = os.path.join(proteome_folder, 'blast_processing')
        ff.create_directory(blast_processing_folder)
        # Create BLAST database files
        [blast_db_files,
        same_protein_other_annotations,
        all_alleles] = create_database_files(proteome_file,
                                            clustering_sim,
                                            clustering_cov,
                                            size_ratio,
                                            blast_processing_folder)
        # Save paths to proteome file paths
        proteomes_data_paths.setdefault(proteome_file_base, [proteome_folder, blast_processing_folder, blast_db_files])
        proteomes_data.append((same_protein_other_annotations, all_alleles))

    [translation_dict,
     reps_ids,
     translations_paths] = sf.translate_schema_loci(schema_directory,
                                                    proteome_matcher_output,
                                                    translation_table,
                                                    run_mode)

    # Import Swiss-Prot and TrEMBL records descriptions
    with open(proteome_files[2], 'rb') as dinfile:
        descriptions: Dict[str, str] = pickle.load(dinfile)

    # For better prints get max length of string
    max_id_length: int = len(max(reps_ids, key=len))

    # Get path to the blastp executable
    blast_exec: str = lf.get_tool_path('blastp')
    self_score_dict: Dict[str, float] = bf.calculate_self_score(translations_paths,
                                                                blast_exec,
                                                                proteome_matcher_output,
                                                                max_id_length,
                                                                cpu)
    merge_files: List[List[str]] = [[], []]
    for i, (file_name, paths) in enumerate(proteomes_data_paths.items()):
        if paths[2] is None:
            prf.print_message(f"Skipping proteome file BLAST: {file_name} due to lack of proteins", "warning")
            continue
        # Get proteome file name without extension
        file_name_without_extension: str = file_name.split('.')[0]
        prf.print_message(f"Processing proteome for: {file_name_without_extension}", "info")
        # Get paths
        [proteome_folder, blast_processing_folder, blast_db_files] = paths
        # Run Blasts and save to file
        best_bsr_values_per_proteome_file: Dict[str, Dict[str, Tuple[str, float]]] = {k: {} for k in proteome_file_ids.keys()}
        best_bsr_values = run_blast_for_proteomes(max_id_length,
                                proteome_file_ids,
                                best_bsr_values_per_proteome_file,
                                blast_processing_folder,
                                translations_paths,
                                blast_db_files,
                                self_score_dict,
                                cpu,
                                bsr)
        
        # Save annotations
        header: str = 'Locus\tUniprot_protein_ID\tUniprot_protein_product\tUniprot_protein_short_name\tUniprot_BSR'
        annotations_file: str = os.path.join(proteome_folder, f"{file_name_without_extension}_annotations.tsv")
        not_matched_or_bsr_failed_loci = set(translations_paths.keys()) - set(best_bsr_values.keys())
        with open(annotations_file, 'w') as at:
            at.write(header + '\n')
            for loci, subject_info in best_bsr_values.items():
                subject_id: str = subject_info[0]
                bsr_value: float = subject_info[1]
                desc: str = descriptions[subject_id]
                lname: str = desc.split(subject_id + ' ')[1].split(' OS=')[0]
                sname: str = desc.split('GN=')[1].split(' PE=')[0]
                # If there is no short name, set it to 'NA'
                if sname == '':
                    sname = 'NA'
                # Write the annotations to the file
                at.write(f"{loci}\t{subject_id}\t{lname}\t{sname}\t{bsr_value}\n")
            # Write loci that did not match or failed the BSR threshold
            for loci in not_matched_or_bsr_failed_loci:
                at.write(f"{loci}\tNA\tNA\tNA\tNA\n")
        # Save annotations file
        merge_files[i].append(annotations_file)

    for i, (proteome_file_id, loci_values) in enumerate(list(best_bsr_values_per_proteome_file.items())):
        (same_protein_other_annotations, all_alleles) = proteomes_data[i]
        if same_protein_other_annotations is None or all_alleles is None:
            continue
        for loci_id, values in list(loci_values.items()):
            # Find all of the proteins that reps representes
            proteinid: str = values[0]
            bsr_value = values[1]
            same_protein: Optional[List[str]] = same_protein_other_annotations.get(proteinid)
            # Check if the protein is being represented by another sequence
            if same_protein:
                # For all of the elements that the representative represents add them to the dict
                for id_ in same_protein:
                    # Get genbank file for that ID
                    proteome_file_id_current: str = itf.identify_string_in_dict_get_key(id_, proteome_file_ids)
                    # If the genbank file is the same as the one that the representative is in, skip (may be copy protein)
                    if proteome_file_id == proteome_file_id_current:
                        continue
                    # Verify if genbank file is in the dict
                    best_bsr_values_per_proteome_file[proteome_file_id_current].setdefault(loci_id, (id_, bsr_value))
            # For clustered elements
            rep_cluster = all_alleles.get(proteinid)
            if rep_cluster:
                for values in rep_cluster:
                    id_ = values[0]
                    # Get genbank file for that ID
                    proteome_file_id_current = itf.identify_string_in_dict_get_key(id_, proteome_file_ids)
                    # If the genbank file is the same as the one that the representative is in, skip (may be paralogous protein)
                    if proteome_file_id == proteome_file_id_current:
                        continue
                    # Verify if genbank file is in the dict
                    best_bsr_values_per_proteome_file[proteome_file_id_current].setdefault(loci_id, (id_, bsr_value))
        
    # Save best annotations per proteome file
    best_annotations_per_proteome_file: str = os.path.join(proteome_matcher_output, "best_annotations_per_proteome_file")
    ff.create_directory(best_annotations_per_proteome_file)
    # Create Swiss-Prot and TrEMBL folders
    swiss_prot_folder: str = os.path.join(best_annotations_per_proteome_file, 'Swiss-Prot')
    ff.create_directory(swiss_prot_folder)
    trembl_folder: str = os.path.join(best_annotations_per_proteome_file, 'TrEMBL')
    ff.create_directory(trembl_folder)

    for file, loci_results in best_bsr_values_per_proteome_file.items():
        # Save what loci each proteome file matched
        matched_loci: Dict[str, List[str]] = {'swiss-prot': [], 'trembl': []}
        # Create Swiss-Prot and TrEMBL annotations files
        swiss_prot_annotations: str = os.path.join(swiss_prot_folder, f"{file}_Swiss-Prot_annotations.tsv")
        trembl_annotations: str = os.path.join(trembl_folder, f"{file}_TrEMBL_annotations.tsv")
        if file in proteome_ids_to_add:
            merge_files[0].append(swiss_prot_annotations)
            merge_files[1].append(trembl_annotations)
        with open(swiss_prot_annotations, 'w') as sp, open(trembl_annotations, 'w') as tr:
            sp.write(header + '\n')
            tr.write(header + '\n')
            for loci, subject_info in loci_results.items():
                subject_id = subject_info[0]
                bsr_value = subject_info[1]
                desc = descriptions[subject_id]
                lname= desc.split(subject_id + ' ')[1].split(' OS=')[0]
                sname = desc.split('GN=')[1].split(' PE=')[0]
                if sname == '':
                    sname = 'NA'
                # Write to the appropriate file based on the start of subject_id
                if subject_id.startswith('sp|'):
                    matched_loci['swiss-prot'].append(loci)
                    sp.write(f"{loci}\t{subject_id}\t{lname}\t{sname}\t{bsr_value}\n")
                elif subject_id.startswith('tr|'):
                    matched_loci['trembl'].append(loci)
                    tr.write(f"{loci}\t{subject_id}\t{lname}\t{sname}\t{bsr_value}\n")
            for proteome_file, loci in matched_loci.items():
                not_matched_or_bsr_failed_loci = set(translations_paths.keys()) - set(loci)   
                for loci in not_matched_or_bsr_failed_loci:
                    if proteome_file == 'swiss-prot':
                        sp.write(f"{loci}\tNA\tNA\tNA\tNA\n")
                    else:
                        tr.write(f"{loci}\tNA\tNA\tNA\tNA\n")

    # Merge all annotations files that user wants
    merged_annotations_file_list = []
    for merge_annotations in merge_files:
        if len(merge_annotations) == 0:
            continue
        merged_annotations_file: str = os.path.join(output_directory, f"best_proteomes_annotations_{'swiss_prot' if i == 0 else 'trEMBL'}.tsv")
        merged_annotations_file_list.append(merged_annotations_file)
        pf.merge_files_into_same_file_by_key(merge_annotations, 'Locus', merged_annotations_file)

    return merged_annotations_file_list
