import os
import pickle
import concurrent.futures
from typing import Dict, List, Set, Tuple, Optional, Union
from itertools import repeat

try:
    from utils import (sequence_functions as sf,
                       file_functions as ff,
                       clustering_functions as cf,
                       blast_functions as bf,
                       linux_functions as lf,
                       alignments_functions as af,
                       iterable_functions as itf)
except ModuleNotFoundError:
    from SchemaRefinery.utils import (sequence_functions as sf,
                                      file_functions as ff,
                                      clustering_functions as cf,
                                      blast_functions as bf,
                                      linux_functions as lf,
                                      alignments_functions as af,
                                      iterable_functions as itf)

def create_database_files(proteome_file: str, clustering_sim: float, clustering_cov: float, size_ratio: float,
                          blast_processing_folder: str) -> Optional[str]:
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
    same_protein_other_annotations: Dict[str, List[Tuple[str, str, str]]] = {}

    print(f"\nExtracting protein from proteome file: {os.path.basename(proteome_file)}...")
    
    # Fetch protein sequences from the proteome file
    fasta_dict: Dict[str, str] = sf.fetch_fasta_dict(proteome_file, False)
    total_proteins: int = len(fasta_dict)
    
    if total_proteins == 0:
        print(f"No proteins found in {os.path.basename(proteome_file)}")
        return None

    # Process each protein sequence
    for id_, sequence in fasta_dict.items():
        protein_hash: str = sf.hash_sequence(sequence)
        if protein_hash in processed_proteins:
            # If the protein has already been processed, save the other annotations that it may have
            same_protein_other_annotations.setdefault(protein_hash, []).append((id_,))
            continue
        else:
            processed_proteins.add(protein_hash)
            all_translation_dict[id_] = sequence

    print(f"\nOut of {total_proteins} protein sequences, {len(all_translation_dict)} are unique proteins\n")
    
    print("Clustering protein sequences...")
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
    
    print(f"Clustered {len(all_translation_dict)} into {len(reps_sequences)} clusters.")
    
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
        
    return blast_db_files

def run_blast_for_proteomes(max_id_length: Dict[str, str], proteome_file_ids: Dict[str, List[str]],
                            best_bsr_values_per_proteome_file: Dict[str, Dict[str, List[Union[str, float]]]],
                            blast_processing_folder: str, translations_paths: Dict[str, str],
                            blast_db_files: str, proteome_folder: str, file_name_without_extension: str,
                            descriptions: Dict[str, str], self_score_dict: Dict[str, float], cpu: int, bsr: float) -> None:
    """
    Run BLAST for proteomes to calculate self-scores and BSR values, and save the annotations.

    Parameters
    ----------
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
    print("Running BLASTp...")
    blastp_results_folder: str = os.path.join(blast_processing_folder, 'blastp_results')
    ff.create_directory(blastp_results_folder)
    
    # Run BLASTp between all BLASTn matches (rep vs all its BLASTn matches).
    bsr_values: Dict[str, Dict[str, float]] = {}
    best_bsr_values: Dict[str, Tuple[str, float]] = {}
    total_blasts: int = len(translations_paths)
    i: int = 1
    
    with concurrent.futures.ProcessPoolExecutor(max_workers=cpu) as executor:
        for res in executor.map(bf.run_blastdb_multiprocessing,
                                repeat(get_blastp_exec),
                                repeat(blast_db_files),
                                translations_paths.values(),
                                translations_paths.keys(),
                                repeat(blastp_results_folder)):
            # Get the alignments
            filtered_alignments_dict: Dict[str, Dict[str, Dict[str, float]]]
            filtered_alignments_dict, _, _, _ = af.get_alignments_dict_from_blast_results(res[1], 0, True, False, True, True, False)

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
                    if current_best_bsr and bsr_value > current_best_bsr[1]:
                        best_bsr_values[loci] = (subject_id, bsr_value)
                    else:
                        best_bsr_values[loci] = (subject_id, bsr_value)
                        
                    # Get best value for genbank file
                    genbank_file: str = itf.identify_string_in_dict_get_key(subject_id, proteome_file_ids)
                    current_best_in_genbank_file: List[Union[str, float]] = best_bsr_values_per_proteome_file[genbank_file].get(loci)
                    if current_best_in_genbank_file and bsr_value > current_best_in_genbank_file[1]:
                        best_bsr_values_per_proteome_file[genbank_file][loci] = [subject_id, bsr_value]
                    else:
                        best_bsr_values_per_proteome_file[genbank_file][loci] = [subject_id, bsr_value]

            print(f"\rRunning BLASTp for cluster representatives matches: {res[0]} - {i}/{total_blasts: <{max_id_length}}", end='', flush=True)
            i += 1

    # Save annotations
    header: str = 'Locus\tProtein_ID\tProtein_product\tProtein_short_name\tBSR'
    annotations_file: str = os.path.join(proteome_folder, f"{file_name_without_extension}_annotations.tsv")
    not_matched_or_bsr_failed_loci = set(translations_paths.keys()) - set(best_bsr_values.keys())
    with open(annotations_file, 'w') as at:
        at.write(header + '\n')
        for loci, subject_info in best_bsr_values.items():
            subject_id: str = subject_info[0]
            split_subject_id: str = subject_id.split('|')[1]
            bsr_value: float = subject_info[1]
            desc: str = descriptions[subject_id]
            lname: str = desc.split(subject_id + ' ')[1].split(' OS=')[0]
            sname: str = desc.split('GN=')[1].split(' PE=')[0]
            # If there is no short name, set it to 'NA'
            if sname == '':
                sname = 'NA'
            # Write the annotations to the file
            at.write(f"{loci}\t{split_subject_id}\t{lname}\t{sname}\t{bsr_value}\n")
        # Write loci that did not match or failed the BSR threshold
        for loci in not_matched_or_bsr_failed_loci:
            at.write(f"{loci}\tNA\tNA\tNA\tNA\n")

    return annotations_file

def proteome_matcher(proteome_files: List[str], proteome_file_ids: Dict[str, List[str]], 
                     schema_directory: str, output_directory: str, cpu: int, bsr: float,
                     translation_table: int, clustering_sim: float,
                     clustering_cov: float, size_ratio: float, run_mode: str) -> None:
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
    proteomes_data_paths: Dict[str, List[str]] = {}
    for proteome_file in proteome_files[:2]:
        # Get proteome file name
        proteome_file_base: str = os.path.basename(proteome_file)
        # Create folder for proteome processing
        proteome_folder: str = os.path.join(proteome_matcher_output, f"{proteome_file_base.split('.')[0]}_processing")
        # Create directory for proteome BLAST processing
        blast_processing_folder: str = os.path.join(proteome_folder, 'blast_processing')
        ff.create_directory(blast_processing_folder)
        # Create BLAST database files
        blast_db_files: Optional[str] = create_database_files(proteome_file,
                                                              clustering_sim,
                                                              clustering_cov,
                                                              size_ratio,
                                                              blast_processing_folder)
        # Save paths to proteome file paths
        proteomes_data_paths.setdefault(proteome_file_base, [proteome_folder, blast_processing_folder, blast_db_files])

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
    
    print('\n')
    # Get path to the blastp executable
    blast_exec: str = lf.get_tool_path('blastp')
    self_score_dict: Dict[str, float] = bf.calculate_self_score(translations_paths,
                                                                blast_exec,
                                                                proteome_matcher_output,
                                                                max_id_length,
                                                                cpu)
    print('\n')
    annotations_files = []
    for file_name, paths in proteomes_data_paths.items():
        if paths[2] is None:
            print(f"\nSkipping proteome file BLAST: {file_name} due to lack of proteins")
            continue

        file_name_without_extension: str = file_name.split('.')[0]
        print(f"\nProcessing proteome for: {file_name_without_extension}")
        # Get paths
        proteome_folder: str
        blast_processing_folder: str
        blast_db_files: str
        [proteome_folder, blast_processing_folder, blast_db_files] = paths
        # Run Blasts and save to file
        best_bsr_values_per_proteome_file: Dict[str, Dict[str, List[Union[str, float]]]] = {k: {} for k in proteome_file_ids.keys()}
        annotations_file = run_blast_for_proteomes(max_id_length,
                                proteome_file_ids,
                                best_bsr_values_per_proteome_file,
                                blast_processing_folder,
                                translations_paths,
                                blast_db_files,
                                proteome_folder,
                                file_name_without_extension,
                                descriptions,
                                self_score_dict,
                                cpu,
                                bsr)
        # Save annotations file
        annotations_files.append(annotations_file)
        
    # Save best annotations per proteome file
    best_annotations_per_proteome_file: str = os.path.join(proteome_matcher_output, "best_annotations_per_proteome_file")
    ff.create_directory(best_annotations_per_proteome_file)
    for file, loci_results in best_bsr_values_per_proteome_file.items():
        annotations_file_proteome: str = os.path.join(best_annotations_per_proteome_file, f"{file}.tsv")
        with open(annotations_file_proteome, 'w') as at:
            for loci, subject_info in loci_results.items():
                subject_id: str = subject_info[0]
                bsr_value: float = subject_info[1]
                desc: str = descriptions[subject_id]
                lname: str = desc.split(subject_id + ' ')[1].split(' OS=')[0]
                sname: str = desc.split('GN=')[1].split(' PE=')[0]
                if sname == '':
                    sname = 'NA'
                # Write the annotations to the file
                at.write(f"{loci}\t{subject_id}\t{lname}\t{sname}\t{bsr_value}\n")

        return annotations_files[0], annotations_files[1]
