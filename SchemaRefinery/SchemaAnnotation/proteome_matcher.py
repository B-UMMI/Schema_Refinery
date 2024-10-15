import os
import pickle
import concurrent.futures
from typing import Dict, List, Set, Tuple, Optional
from itertools import repeat

try:
    from utils import (sequence_functions as sf,
                       file_functions as ff,
                       clustering_functions as cf,
                       blast_functions as bf,
                       linux_functions as lf,
                       alignments_functions as af)
except ModuleNotFoundError:
    from SchemaRefinery.utils import (sequence_functions as sf,
                                      file_functions as ff,
                                      clustering_functions as cf,
                                      blast_functions as bf,
                                      linux_functions as lf,
                                      alignments_functions as af)

def create_database_files(proteome_file: str, clustering_sim: float, clustering_cov: float, size_ratio: float,
                          output_directory: str) -> Optional[str]:
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
    output_directory : str
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
    
    print(f"Clustered {len(all_translation_dict)} into {len(reps_sequences)} clusters.\n")
    
    # Save clustered protein sequences to a file
    clustered_protein_master_file: str = os.path.join(output_directory, f"{file_name}.fasta")
    with open(clustered_protein_master_file, 'w') as outfile:
        for protein_id, values in reps_sequences.items():
            outfile.write(f">{protein_id}\n{values}\n")
    
    # Create BLAST database
    blastdb_path: str = os.path.join(output_directory, 'blastdb')
    ff.create_directory(blastdb_path)
    blast_db_files: str = os.path.join(blastdb_path, 'genbank_protein_db')
    makeblastdb_exec: str = lf.get_tool_path('makeblastdb')
    bf.make_blast_db(makeblastdb_exec, clustered_protein_master_file, blast_db_files, 'prot')
        
    return blast_db_files

def self_score_calculation(reps_ids: Dict[str, str], self_score_folder: str,
                           translations_paths: Dict[str, str], cpu: int) -> Dict[str, float]:
    """
    Calculate self-scores for representative protein sequences using BLASTp.

    Parameters
    ----------
    reps_ids : Dict[str, str]
        Dictionary of representative IDs.
    self_score_folder : str
        Path to the folder for storing self-score results.
    translations_paths : Dict[str, str]
        Dictionary of translation paths.
    cpu : int
        Number of CPUs to use.

    Returns
    -------
    Dict[str, float]
        Dictionary of self-scores for each representative ID.
    """
    # Create directory for self-score results
    print("\nCalculate self-score for the schema...")
    self_score_reps_folder: str = os.path.join(self_score_folder, 'self_score_reps')
    ff.create_directory(self_score_reps_folder)
    
    # Initialize dictionary to store self-scores
    self_score_dict: Dict[str, float] = {}
    max_id_length: int = len(max(reps_ids, key=len))
    
    # Get path to the blastp executable
    get_blastp_exec: str = lf.get_tool_path('blastp')
    i: int = 1
    self_score: float

    # Calculate self-scores using multiprocessing
    with concurrent.futures.ProcessPoolExecutor(max_workers=cpu) as executor:
        for res in executor.map(bf.run_self_score_multiprocessing,
                                translations_paths.keys(),
                                repeat(get_blastp_exec),
                                translations_paths.values(),
                                repeat(self_score_reps_folder)):
            
            # Extract self-scores from BLAST results
            _, self_score, _, _ = af.get_alignments_dict_from_blast_results(res[1], 0, True, True, True, True, True)
    
            # Update self-score dictionary
            self_score_dict.update(self_score)
                            
            print(f"\rRunning BLASTp to calculate self-score for {res[0]: <{max_id_length}}", end='', flush=True)
            i += 1

    return self_score_dict

def run_blast_for_proteomes(reps_ids: Dict[str, str], blast_processing_folder: str, translations_paths: Dict[str, str],
                            blast_db_files: str, proteome_folder: str, file_name_without_extension: str,
                            descriptions: Dict[str, str], self_score_dict: Dict[str, float], cpu: int, bsr: float) -> None:
    """
    Run BLAST for proteomes to calculate self-scores and BSR values, and save the annotations.

    Parameters
    ----------
    reps_ids : Dict[str, str]
        Dictionary of representative IDs.
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
    
    # Get the maximum length of the IDs for better prints
    max_id_length: int = len(max(reps_ids, key=len))
    
    # Run BLASTp
    print("\nRunning BLASTp...")
    blastp_results_folder: str = os.path.join(blast_processing_folder, 'blastp_results')
    ff.create_directory(blastp_results_folder)
    
    # Run BLASTp between all BLASTn matches (rep vs all its BLASTn matches).
    bsr_values: Dict[str, Dict[str, float]] = {}
    best_bsr_values: Dict[str, Tuple[str, float]] = {}
    total_blasts: int = len(reps_ids)
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

            print(f"\rRunning BLASTp for cluster representatives matches: {res[0]} - {i}/{total_blasts: <{max_id_length}}", end='', flush=True)
            i += 1

    # Save annotations
    header: str = 'Locus\tProtein_ID\tProtein_product\tProtein_short_name\tBSR'
    annotations_file: str = os.path.join(proteome_folder, f"{file_name_without_extension}_annotations.tsv")
    with open(annotations_file, 'w') as at:
        at.write(header + '\n')
        for loci, subject_info in best_bsr_values.items():
            subject_id: str = subject_info[0]
            split_subject_id: str = subject_id.split('|')[1]
            bsr_value: float = subject_info[1]
            desc: str = descriptions[subject_id]
            lname: str = desc.split(subject_id + ' ')[1].split(' OS=')[0]
            sname: str = desc.split('GN=')[1].split(' PE=')[0]
            # Write the annotations to the file
            at.write(f"{loci}\t{split_subject_id}\t{lname}\t{sname}\t{bsr_value}\n")

def proteome_matcher(proteome_files: List[str], schema_directory: str,
                     output_directory: str, cpu: int, bsr: float,
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
    # Create BLAST database files for each proteome
    proteomes_data_paths: Dict[str, List[str]] = {}
    for proteome_file in proteome_files[:2]:
        # Get proteome file name
        proteome_file_base: str = os.path.basename(proteome_file)
        # Create folder for proteome processing
        proteome_folder: str = os.path.join(output_directory, f"{proteome_file_base.split('.')[0]}_processing")
        # Create directory for proteome BLAST processing
        blast_processing_folder: str = os.path.join(proteome_folder, 'blast_processing')
        ff.create_directory(blast_processing_folder)
        # Create BLAST database files
        blast_db_files: Optional[str] = create_database_files(proteome_file, clustering_sim, clustering_cov, size_ratio, blast_processing_folder)
        # Save paths to proteome file paths
        proteomes_data_paths.setdefault(proteome_file, [proteome_folder, blast_processing_folder, blast_db_files])

    # BLAST alleles for each locus against file with all CDSs from origin genomes
    fasta_files_dict: Dict[str, str] = {
        loci.split('.')[0]: os.path.join(schema_directory, loci)
        for loci in os.listdir(schema_directory)
        if os.path.isfile(os.path.join(schema_directory, loci)) and loci.endswith('.fasta')
    }
    # Short folder
    short_folder: str = os.path.join(schema_directory, 'short')
    fasta_files_short_dict: Dict[str, str] = {
        loci.split('.')[0].split('_')[0]: os.path.join(short_folder, loci)
        for loci in os.listdir(short_folder)
        if os.path.isfile(os.path.join(short_folder, loci)) and loci.endswith('.fasta')
    }
    
    # Translate sequences and save to FASTA file
    files_to_run: Dict[str, str] = fasta_files_dict if run_mode == 'alleles' else fasta_files_short_dict
    print("Translating sequences...")
    reps_translations_folder: str = os.path.join(output_directory, 'reps_translations')
    ff.create_directory(reps_translations_folder)
    translation_dict: Dict[str, str] = {}
    reps_ids: Dict[str, List[str]] = {}
    translations_paths: Dict[str, str] = {}
    for loci, loci_path in files_to_run.items():
        fasta_dict: Dict[str, str] = sf.fetch_fasta_dict(loci_path, False)
        print(f"\rTranslating schema reps: {loci}", end='', flush=True)
        for allele_id, sequence in fasta_dict.items():
            reps_ids.setdefault(loci, []).append(allele_id)
            
            # Translate sequences and update translation dictionary
            trans_path_file: str = os.path.join(reps_translations_folder, f"{loci}.fasta")
            translations_paths[loci] = trans_path_file
            trans_dict: Dict[str, str]
            trans_dict, _, _ = sf.translate_seq_deduplicate(fasta_dict,
                                                            trans_path_file,
                                                            None,
                                                            0,
                                                            False,
                                                            translation_table,
                                                            False)
            translation_dict.update(trans_dict)

    # Import Swiss-Prot and TrEMBL records descriptions
    with open(proteome_files[2], 'rb') as dinfile:
        descriptions: Dict[str, str] = pickle.load(dinfile)

    # Calculate self-score for representatives
    self_score_folder: str = os.path.join(output_directory, 'self_score')
    self_score_dict: Dict[str, float] = self_score_calculation(reps_ids, self_score_folder, translations_paths, cpu)

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
        run_blast_for_proteomes(reps_ids,
                                blast_processing_folder,
                                translations_paths,
                                blast_db_files,
                                proteome_folder,
                                file_name_without_extension,
                                descriptions,
                                self_score_dict,
                                cpu,
                                bsr)