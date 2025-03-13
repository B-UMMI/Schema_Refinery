import os
import re
from Bio import SeqIO
from typing import Dict, List, Tuple

try:
    from utils import (
                        sequence_functions as sf,
                        blast_functions as bf,
                        linux_functions as lf,
                        file_functions as ff,
                        alignments_functions as af,
                        Types as tp,
                        print_functions as pf,
                        logger_functions as logf,
                        globals as gb
    )
except ModuleNotFoundError:
    from SchemaRefinery.utils import (
                                    sequence_functions as sf,
                                    blast_functions as bf,
                                    linux_functions as lf,
                                    file_functions as ff,
                                    alignments_functions as af,
                                    Types as tp,
                                    print_functions as pf,
                                    logger_functions as logf,
                                    globals as gb
                                    
    )

def get_schema_files(schema_directory: str) -> Tuple[Dict[str, str], Dict[str, str]]:
    """
    Identify all of the FASTA files in the schema directory and its 'short' subdirectory.

    Parameters
    ----------
    schema_directory : str
        Path to the directory containing schema FASTA files.

    Returns
    -------
    Tuple[Dict[str, str], Dict[str, str]]
        A tuple containing two dictionaries:
        - The first dictionary maps loci names to their file paths in the schema directory.
        - The second dictionary maps loci names to their file paths in the 'short' subdirectory.
    """
    # Identify all of the FASTA files in the schema directory
    fasta_files_dict: Dict[str, str] = {
        loci.split('.')[0]: os.path.join(schema_directory, loci)
        for loci in os.listdir(schema_directory)
        if os.path.isfile(os.path.join(schema_directory, loci)) and loci.endswith('.fasta')
    }
    
    # Identify all of the FASTA files in the 'short' subdirectory
    short_folder: str = os.path.join(schema_directory, 'short')
    fasta_files_short_dict: Dict[str, str] = {
        loci.split('.')[0].split('_')[0]: os.path.join(short_folder, loci)
        for loci in os.listdir(short_folder)
        if os.path.isfile(os.path.join(short_folder, loci)) and loci.endswith('.fasta')
    }
    
    return fasta_files_dict, fasta_files_short_dict


def run_blasts_match_schemas(query_translations_paths: Dict[str, str], blast_db_files: str,
                             blast_folder: str, self_score_dict: Dict[str, float], max_id_length: int,
                             get_blastp_exec: str, bsr: float, cpu: int) -> Dict[str, Dict[str, float]]:
    """
    Run BLASTp to match schemas and calculate BSR values.

    Parameters
    ----------
    query_translations_paths : dict
        Dictionary with keys as query identifiers and values as paths to the query translation files.
    blast_db_files : str
        Path to the BLAST database files.
    blast_folder : str
        Path to the folder where BLAST results will be stored.
    self_score_dict : dict
        Dictionary with query identifiers as keys and their self-scores as values.
    max_id_length : int
        Maximum length of the query identifiers.
    get_blastp_exec : str
        Path to the BLASTp executable.
    bsr : float
        BSR threshold value.
    cpu : int
        Number of CPU cores to use for multiprocessing.

    Returns
    -------
    dict
        Dictionary with loci identifiers as keys and tuples of the best subject identifier and BSR value as values.
    """
    
    # Run BLASTp
    pf.print_message("Running BLASTp...", "info")
    blastp_results_folder: str = os.path.join(blast_folder, 'blastp_results')
    ff.create_directory(blastp_results_folder)
    
    # Initialize dictionaries to store BSR values
    bsr_values: Dict[str, Dict[str, float]] = {}
    best_bsr_values: Dict[str, Dict[str, float]] = {}
    total_blasts: int = len(query_translations_paths)
    blastp_results_files = bf.run_blastp_operations(cpu,
                                                    get_blastp_exec,
                                                    blast_db_files,
                                                    query_translations_paths,
                                                    blastp_results_folder,
                                                    total_blasts,
                                                    max_id_length)

    for blast_result_file in blastp_results_files:
        # Get the alignments
        filtered_alignments_dict: tp.BlastDict 
        filtered_alignments_dict, _, _ = af.get_alignments_dict_from_blast_results_simplified(blast_result_file,
                                                                                                    0,
                                                                                                    False,
                                                                                                    True)
       

        # Since BLAST may find several local alignments, choose the largest one to calculate BSR.
        for query, subjects_dict in filtered_alignments_dict.items():
            # Get the loci name
            query_loci_id: str = query.split('_')[0]

            best_bsr_values.setdefault(query_loci_id, {})
            # Create the dict of the query
            bsr_values.setdefault(query, {})
            for subject_id, results in subjects_dict.items():
                subject_loci_id = subject_id.split('_')[0]
                # Highest score (First one)
                subject_score: float = next(iter(results.values()))['score']
                # Calculate BSR value
                computed_score: float = bf.compute_bsr(subject_score, self_score_dict[query])
                # Check if the BSR value is higher than the threshold
                if computed_score >= bsr:
                    # Round BSR values if they are superior to 1.0 to 1 decimal place
                    if computed_score > 1.0:
                        computed_score = round(computed_score, 1)
                    # Save all of the different matches that this query had and their BSR values
                    bsr_values[query].update({subject_id: computed_score})
                else:
                    continue
                # Save the best match for each query and subject matches
                subject_loci_id = subject_id.split('_')[0]
                if not best_bsr_values[query_loci_id].get(subject_loci_id):
                    best_bsr_values[query_loci_id][subject_loci_id] = computed_score
                elif computed_score > best_bsr_values[query_loci_id][subject_loci_id]:
                    best_bsr_values[query_loci_id][subject_loci_id] = computed_score

    best_bsr_values = {query: match for query, match in best_bsr_values.items() if match}

    return best_bsr_values


def write_best_blast_matches_to_file(best_bsr_values: Dict[str, Dict[str, float]],
                                     query_translations_paths: Dict[str, str], 
                                     subject_translations_paths: Dict[str, str], 
                                     output_folder: str, 
                                     rep_vs_alleles: bool, 
                                     process_name: str) -> str:
    """
    Write the best BLAST matches to a file.

    Parameters
    ----------
    best_bsr_values : dict
        Dictionary with loci identifiers as keys and tuples of the best subject identifier and BSR value as values.
    query_translations_paths : dict
        Dictionary with keys as query identifiers and values as paths to the query translation files.
    subject_translations_paths : dict
        Dictionary with keys as subject identifiers and values as paths to the subject translation files.
    output_folder : str
        Path to the folder where the output file will be stored.

    Returns
    -------
    None
    """

    # Path to output files
    best_blast_matches_file = os.path.join(output_folder, "Match_Schemas_Results.tsv")
    existing_matches_file = os.path.join(output_folder, "existing_matches.txt")

    # Load existing matches from existing ,atches file to avoid repetition across runs
    existing_matches = set()
    written_queries = set()

    if os.path.exists(existing_matches_file):
        with open(existing_matches_file, "r") as f:
            for line in f:
                parts = line.strip().split("\t")
                if parts:  
                    existing_matches.add(line.strip())
                    written_queries.add(parts[0])

    # Initialize lists
    matched_entries = []
    non_matched_query = []
    non_matched_subject = []

    # Identify loci that were not matched
    not_matched_loci = [query for query in query_translations_paths.keys() if query not in best_bsr_values.keys()]
    not_matched_subject = [subject for subject in subject_translations_paths.keys() if subject not in {s for d in best_bsr_values.values() for s in d.keys()}]

    # Check if the best matches file exists
    file_exists = os.path.exists(best_blast_matches_file)

    # Read existing file to avoid duplicates in this run
    if file_exists:
        with open(best_blast_matches_file, "r") as f:
            next(f)
            for line in f:
                parts = line.strip().split("\t")
                if parts:
                    existing_matches.add(line.strip())
                    written_queries.add(parts[0])

    # Write best matches
    with open(best_blast_matches_file, "a" if file_exists else "w") as out:
        if not file_exists:
            out.write("Query\tSubject\tBSR\tProcess\n")

        for query, match in best_bsr_values.items():
            if query in written_queries:
                continue  # Skip if query was already written
            for subject, computed_score in match.items():
                entry = f"{query}\t{subject}\t{computed_score}"
                if entry not in existing_matches:
                    existing_matches.add(f"{entry}")
                    matched_entries.append((query, f"{entry}\t{process_name}"))
            written_queries.add(query)

        # Process unmatched queries
        for query in not_matched_loci:
            if query in written_queries:
                continue  # Skip if query was already written
            entry = f"{query}\tNot matched\tNA"
            non_matched_query.append((query, f"{entry}\t{process_name}"))
            written_queries.add(query) 

        # Process unmatched subjects
        for subject in not_matched_subject:
            entry = f"Not matched\t{subject}\tNA"
            non_matched_subject.append((subject, f"{entry}\t{process_name}"))

        # Sort and write all entries
        for _, entry in sorted(matched_entries, key=lambda x: x[0]):
            out.write(entry + "\n")
        if rep_vs_alleles:
            for _, entry in sorted(non_matched_query, key=lambda x: x[0]):
                if process_name == "rep_vs_alleles":
                    out.write(entry + "\n")
            for _, entry in sorted(non_matched_subject, key=lambda x: x[0]):
                if process_name == "rep_vs_alleles":
                    out.write(entry + "\n")
        else:
            for _, entry in sorted(non_matched_query, key=lambda x: x[0]):
                if process_name == "rep_vs_rep":
                    out.write(entry + "\n")
            for _, entry in sorted(non_matched_subject, key=lambda x: x[0]):
                if process_name == "rep_vs_rep":
                    out.write(entry + "\n")
    

    # Save `existing_matches` to a file for persistence
    with open(existing_matches_file, "w") as f:
        for entry in existing_matches:
            f.write(entry + "\n")

    return best_blast_matches_file



def match_schemas(first_schema_directory: str, second_schema_directory: str, output_directory: str, bsr: float,
                  translation_table: int, cpu: int, no_cleanup: bool, rep_vs_alleles: bool) -> str:
    """
    Match schemas between query and subject directories.

    Parameters
    ----------
    first_schema_directory : str
        Path to the first schema directory.
    second_schema_directory : str
        Path to the second schema directory.
    output_directory : str
        Path to the output directory.
    bsr : float
        BLAST Score Ratio value.
    translation_table : int
        Genetic code used for translation.
    cpu : int
        Number of CPU cores to use.
    no_cleanup : bool
        If True, temporary files will not be removed.
    rep_vs_alleles: bool
        If True then after the rep vs rep Blast the program will run a second Blast with rep vs alleles.

    Returns
    -------
    None
    """
    # A schema files
    a_files: Dict[str, str]
    a_files_short: Dict[str, str]
    a_files, a_files_short = get_schema_files(first_schema_directory)
    # B schema files
    b_files: Dict[str, str]
    b_files_short: Dict[str, str]
    b_files, b_files_short = get_schema_files(second_schema_directory)

    # Choose which schema will be the query (the one with the higher average of alleles per loci)
    total_alleles_a = 0
    avg_a = 0

    for query_loci, fasta_path in a_files.items():
        allele_dict: Dict[str, str] = sf.fetch_fasta_dict(fasta_path, False)
        num_alleles = len(allele_dict)
        total_alleles_a += num_alleles
        avg_a = total_alleles_a/(len(a_files))

    pf.print_message(f"Total alleles in First Schema: {total_alleles_a}. Total Loci: {len(a_files)}. And an average of {round(avg_a)} alleles per loci.", "info")


    total_alleles_b = 0
    avg_b = 0
    for query_loci, fasta_path in b_files.items():
        allele_dict: Dict[str, str] = sf.fetch_fasta_dict(fasta_path, False)
        num_alleles = len(allele_dict)
        total_alleles_b += num_alleles
        avg_b = total_alleles_b/(len(b_files))

    pf.print_message(f"Total alleles in Second Schema: {total_alleles_b}. Total Loci: {len(b_files)}. And an average of {round(avg_b)} alleles per loci.", "info")


    # Set first or second schema as the query or subject
    if avg_a >= avg_b:
        query_fastas_hash: Dict[str, str] = a_files
        subject_fastas_hash: Dict[str, str] = b_files
        query_fastas_rep: Dict[str, str] = a_files_short
        subject_fastas_rep: Dict[str, str] = b_files_short
        pf.print_message(f"{first_schema_directory} set as Query.", "info")
        pf.print_message(f"{second_schema_directory} set as Subject.", "info")
        pf.print_message(f"Total alleles in Query Schema: {total_alleles_a}. Total Loci: {len(a_files)}. And an average of {round(avg_a)} alleles per loci.", "info")
        pf.print_message(f"Total alleles in Subject Schema: {total_alleles_b}. Total Loci: {len(b_files)}. And an average of {round(avg_b)} alleles per loci.", "info")

    else:
        query_fastas_hash: Dict[str, str] = b_files
        subject_fastas_hash: Dict[str, str] = a_files
        query_fastas_rep: Dict[str, str] = b_files_short
        subject_fastas_rep: Dict[str, str] = a_files_short
        pf.print_message(f"{second_schema_directory} set as Query.", "info")
        pf.print_message(f"{first_schema_directory} set as Subject.", "info")
        pf.print_message(f"Total alleles in Query Schema: {total_alleles_b}. Total Loci: {len(b_files)}. And an average of {round(avg_b)} alleles per loci.", "info")
        pf.print_message(f"Total alleles in Subject Schema: {total_alleles_a}. Total Loci: {len(a_files)}. And an average of {round(avg_a)} alleles per loci.", "info")

    # Create the output directories
    blast_folder: str = os.path.join(output_directory, 'blast_processing')
    ff.create_directory(blast_folder)
    # Directories for complete schemas
    query_translation_folder: str = os.path.join(blast_folder, 'Query_Translation')
    ff.create_directory(query_translation_folder)
    subject_translation_folder: str = os.path.join(blast_folder, 'Subject_Translation')
    ff.create_directory(subject_translation_folder)
    query_untranslation_folder: str = os.path.join(blast_folder, 'Query_Not_Translated')
    ff.create_directory(query_untranslation_folder)
    # Directories for the representatives
    query_translation_rep_folder: str = os.path.join(blast_folder, 'Query_Translation_Rep')
    ff.create_directory(query_translation_rep_folder)
    subject_translation_rep_folder: str = os.path.join(blast_folder, 'Subject_Translation_Rep')
    ff.create_directory(subject_translation_rep_folder)
    query_untranslation_rep_folder: str = os.path.join(blast_folder, 'Query_Not_Translated_Rep')
    ff.create_directory(query_untranslation_rep_folder)

    deleted_files = 0 

    len_query_fastas: int = len(query_fastas_hash)
    len_subject_fasta: int = len(subject_fastas_hash)
    len_query_rep_fastas: int = len(query_fastas_rep)
    len_subject_rep_fasta: int = len(subject_fastas_rep)

    # Process query FASTA files for representatives
    query_translation_rep_dict: Dict[str, str] = {}
    query_ids: Dict[str, List[List[str]]] = {}
    query_translations_rep_paths: Dict[str, str] = {}
    i = 0
    pf.print_message("", 'info')
    pf.print_message("Translating sequences for representative query schema...", "info")
    for query_loci, path in query_fastas_rep.items():
        i += 1
        pf.print_message(f"Translated representative query loci FASTA: {i}/{len_query_rep_fastas}", "info", end='\r', flush=True)
        # Get the fasta sequences for the query
        fasta_dict: Dict[str, str] = sf.fetch_fasta_dict(path, False)
        # Save the IDs of the alleles
        query_ids.setdefault(query_loci, []).append([allele_id for allele_id in fasta_dict.keys()])
        # Create translation file path
        query_fasta_rep_translation = os.path.join(query_translation_rep_folder, f"{query_loci}_rep_translation.fasta")
        query_untranslated_rep_seq = os.path.join(query_untranslation_rep_folder, f"{query_loci}_rep_not_translated.fasta")
        # Translate sequences, save untranslated sequences and update translation dictionary
        query_translations_rep_paths[query_loci] = query_fasta_rep_translation
        trans_dict, _, untras_seq= sf.translate_seq_deduplicate(fasta_dict,
                                                        query_fasta_rep_translation,
                                                        query_untranslated_rep_seq,
                                                        0,
                                                        False,
                                                        translation_table,
                                                        False)
        # Update the query translation and protein dictionaries
        query_translation_rep_dict.update(trans_dict)

        # Delete the files that are empty --> the translation has failed
        for file_name in os.listdir(query_translation_rep_folder):
            file_path = os.path.join(query_translation_rep_folder, file_name)
            if os.path.isfile(file_path) and os.path.getsize(file_path) == 0:
                os.remove(file_path)
                seq_id = file_name.replace("_rep_translation.fasta", "")
                query_translations_rep_paths.pop(seq_id, None)
                deleted_files += 1
                pf.print_message(f"Deleted empty file: {file_path}", "info", end='\r', flush=True)
    # If all the sequences are translated the folders are deleted
    if os.path.isdir(query_untranslation_rep_folder) and not os.listdir(query_untranslation_rep_folder):
        os.rmdir(query_untranslation_rep_folder)
    pf.print_message(f"Cleanup complete: Removed {deleted_files} empty translation files.", "info")


    # Process subject FASTA files for representatives
    subject_translation_rep_dict: Dict[str, str] = {}
    subject_ids: Dict[str, List[List[str]]] = {}
    subject_translations_rep_paths: Dict[str, str] = {}
    master_file_rep_path: str = os.path.join(blast_folder, 'subject_master_rep_file.fasta') #file with the translation of the representatives
    i = 0
    pf.print_message("Translating sequences for representative subject schema...", "info")
    for subject_loci, path in subject_fastas_rep.items():
        i += 1
        pf.print_message(f"Translated representative subject loci FASTA: {i}/{len_subject_rep_fasta}", "info", end='\r', flush=True)
        # Get the fasta sequences for the subject
        fasta_dict = sf.fetch_fasta_dict(path, False)
        # Save the IDs of the alleles
        subject_ids.setdefault(subject_loci, []).append([allele_id for allele_id in fasta_dict.keys()])
        # Create translation file path
        subject_fasta_rep_translation = os.path.join(subject_translation_rep_folder, f"{subject_loci}_rep_translation.fasta")
        # Translate sequences and update translation dictionary
        subject_translations_rep_paths[subject_loci] = subject_fasta_rep_translation
        trans_dict, _, _= sf.translate_seq_deduplicate(fasta_dict,
                                                        subject_fasta_rep_translation,
                                                        None,
                                                        0,
                                                        False,
                                                        translation_table,
                                                        False)

        # Update the subject translation and proteins dictionaries
        subject_translation_rep_dict.update(trans_dict)


    # Process query FASTA files for the complete query schema
    query_translation_dict: Dict[str, str] = {}
    query_ids: Dict[str, List[List[str]]] = {}
    query_translations_paths: Dict[str, str] = {}
    query_prot_hash: Dict[str, str] = {}
    i = 0
    pf.print_message("", "info")
    pf.print_message("Translating sequences for complete query schema...", "info")
    for query_loci, path in query_fastas_hash.items():
        i += 1
        pf.print_message(f"Translated complete query loci FASTA: {i}/{len_query_fastas}", "info", end='\r', flush=True)
        # Get the fasta sequences for the query
        fasta_dict: Dict[str, str] = sf.fetch_fasta_dict(path, False)
        # Save the IDs of the alleles
        query_ids.setdefault(query_loci, []).append([allele_id for allele_id in fasta_dict.keys()])
        # Create translation file path
        query_fasta_translation = os.path.join(query_translation_folder, f"{query_loci}_translation.fasta")
        query_untranslated_seq = os.path.join(query_untranslation_folder, f"{query_loci}_not_translated.fasta")
        # Translate sequences, save untranslated sequences, get the protein hashes and update translation dictionary
        query_translations_paths[query_loci] = query_fasta_translation
        trans_dict, prot_hashes, untras_seq= sf.translate_seq_deduplicate(fasta_dict,
                                                        query_fasta_translation,
                                                        query_untranslated_seq,
                                                        0,
                                                        False,
                                                        translation_table,
                                                        True)
        # Update the query translation and protein dictionaries
        query_translation_dict.update(trans_dict)
        query_prot_hash.update(prot_hashes)

        # Delete the files that are empty --> the translation has failed
        for file_name in os.listdir(query_translation_folder):
            file_path = os.path.join(query_translation_folder, file_name)
            if os.path.isfile(file_path) and os.path.getsize(file_path) == 0:
                os.remove(file_path)
                seq_id = file_name.replace("_translation.fasta", "")
                query_translations_paths.pop(seq_id, None)
                
                deleted_files += 1
                pf.print_message(f"Deleted empty file: {file_path}", "info", end='\r', flush=True)
    # Delete the folder if all proteins are translated    
    if os.path.isdir(query_untranslation_folder) and not os.listdir(query_untranslation_folder):
        os.rmdir(query_untranslation_folder)
        
    pf.print_message(f"Cleanup complete: Removed {deleted_files} empty translation files.", "info")


    # Process subject FASTA files for the complete subject schema
    subject_translation_dict: Dict[str, str] = {}
    subject_ids: Dict[str, List[List[str]]] = {}
    subject_translations_paths: Dict[str, str] = {}
    subject_prot_hash: Dict[str, str] = {}
    master_file_path: str = os.path.join(blast_folder, 'subject_master_file.fasta') #file with all the translated sequences of the subject schema
    i = 0
    pf.print_message("Translating sequences for complete subject schema...", "info")
    for subject_loci, path in subject_fastas_hash.items():
        i += 1
        pf.print_message(f"Translated complete subject loci FASTA: {i}/{len_subject_fasta}", "info", end='\r', flush=True)
        # Get the fasta sequences for the subject
        fasta_dict = sf.fetch_fasta_dict(path, False)
        # Save the IDs of the alleles
        subject_ids.setdefault(subject_loci, []).append([allele_id for allele_id in fasta_dict.keys()])
        # Create translation file path
        subject_fasta_translation = os.path.join(subject_translation_folder, f"{subject_loci}_translation.fasta")
        # Translate sequences, get protein hashes and update translation dictionary
        subject_translations_paths[subject_loci] = subject_fasta_translation
        trans_dict, prot_hash, _= sf.translate_seq_deduplicate(fasta_dict,
                                                        subject_fasta_translation,
                                                        None,
                                                        0,
                                                        False,
                                                        translation_table,
                                                        True)

        # Update the subject translation and proteins dictionaries
        subject_translation_dict.update(trans_dict)
        subject_prot_hash.update(prot_hash)


    # -------------------------------------------------------------------
    # Comparision of the Query and Subject hashes (the BSR = 1.0)
    # -------------------------------------------------------------------
    # Prepare best BSR values and query translations

    pf.print_message("", "info")
    pf.print_message("Matching hashes between query and subject schema", "info")
    pf.print_message(f"The query schema has {len(query_prot_hash)} hashes.", "info")
    pf.print_message(f"The subject schema has {len(subject_prot_hash)} hashes.", "info")
    best_bsr_values = {}

    # Find common keys (matching protein hashes)
    common_keys = set(query_prot_hash.keys()) & set(subject_prot_hash.keys())
    seen_pairs = set()
    subject_base_dict = {}
    locus_removal = 0

    for prot_hash in common_keys:
        query_ids = query_prot_hash[prot_hash]
        subject_ids = subject_prot_hash[prot_hash]
        
        for query_id in query_ids:
            # Go from allele to loci
            query_base = re.sub(r'_\*?\d+$', '', query_id)
            for subject_id in subject_ids:
                subject_base = re.sub(r'_\*?\d+$', '', subject_id)
                # Only add unique pairs
                if (query_base, subject_base) not in seen_pairs:
                    best_bsr_values.setdefault(query_base, {})[subject_base] = 1.0
                    seen_pairs.add((query_base, subject_base))
                    subject_translations_paths.pop(subject_base, None)
                    subject_translations_rep_paths.pop(subject_base, None)
                    locus_removal += 1

    # Write results to the best matches file
    pf.print_message("Writting results to the output file...", "info")
    best_blast_matches_file = write_best_blast_matches_to_file(best_bsr_values, query_translations_rep_paths, subject_translations_paths, output_directory, rep_vs_alleles,'hashes_vs_hashes')

    # Print out stats
    pf.print_message(f"From the hash comparision {len(seen_pairs)} matches were found.", "info")
    pf.print_message(f"{locus_removal} loci were removed.", "info")
    pf.print_message(f"{len(subject_translations_rep_paths)} representative subject loci will pass for the next match analysis.", "info")

    # Write the master file without the subject sequences that have already be matched
    pf.print_message("Writting rep master file...", "info")
    with open(master_file_rep_path, 'w') as master:
        for locus_id, fasta_path in subject_translations_rep_paths.items():
            with open(fasta_path, 'r') as fasta_file:
                for record in SeqIO.parse(fasta_file, 'fasta'):
                    master.write(f">{record.id}\n{record.seq}\n")


    # -------------------------------------------------------------------
    # Blast with rep vs rep
    # -------------------------------------------------------------------
    # Get Path to the blastp executable
    pf.print_message("", 'info')
    get_blastp_exec: str = lf.get_tool_path('blastp')
    
    # Get the maximum length of the IDs for better prints
    max_id_length: int = len(max(query_translations_rep_paths.keys(), key=len))

    # Calculate self-scores for each query
    self_score_dict: Dict[str, float] = bf.calculate_self_score(query_translations_rep_paths,
                                                              get_blastp_exec,
                                                              blast_folder,
                                                              max_id_length,
                                                              cpu)

    # Create BLAST database
    pf.print_message("Creating Blast database with representatives from subject schema...", "info")
    blastdb_path: str = os.path.join(blast_folder, 'blastdb')
    ff.create_directory(blastdb_path)
    blast_db_files: str = os.path.join(blastdb_path, 'genbank_protein_db')
    makeblastdb_exec: str = lf.get_tool_path('makeblastdb')
    bf.make_blast_db(makeblastdb_exec, master_file_rep_path, blast_db_files, 'prot')

    # Run BLAST rep_vs_rep
    pf.print_message("Running BLASTs (rep vs rep) between schemas...", "info")
    best_bsr_values: Dict[str, Dict[str, float]] = run_blasts_match_schemas(query_translations_rep_paths,
                                                                            blast_db_files,
                                                                            blast_folder,
                                                                            self_score_dict,
                                                                            max_id_length,
                                                                            get_blastp_exec,
                                                                            bsr,
                                                                            cpu)
    pf.print_message("", None)


    # Write the best BLAST matches to a file
    pf.print_message("Writting results to the output file...", "info")
    best_blast_matches_file = write_best_blast_matches_to_file(best_bsr_values, query_translations_rep_paths, subject_translations_rep_paths, output_directory, rep_vs_alleles, 'rep_vs_rep')

    # Remove matched id from the full subject schema file
    locus_removal = 0
    subject_base_list = set()
    for query, match in best_bsr_values.items():
        for subject, score in match.items():
            subject_base = re.sub(r'_\*?\d+$', '', subject)
            subject_base_list.add(subject_base)
    
    for subject in subject_base_list:
        subject_translations_paths.pop(subject, None)
        locus_removal += 1
       
    # Write the sequences to the full master file
    pf.print_message("Writting master file...", "info")
    with open(master_file_path, 'w') as master:
        for locus_id, fasta_path in subject_translations_paths.items():
            with open(fasta_path, 'r') as fasta_file:
                for record in SeqIO.parse(fasta_file, 'fasta'):
                    master.write(f">{record.id}\n{record.seq}\n")

    pf.print_message(f"From the rep vs rep Blast {len(best_bsr_values)} matches were found.", "info")
    pf.print_message(f"{locus_removal} loci were removed.", "info")
    pf.print_message(f"{len(subject_translations_paths)} subject loci have not found a match.", "info")

    # -------------------------------------------------------------------
    # Blast with rep vs alleles
    # -------------------------------------------------------------------
    # Create BLAST database
    if rep_vs_alleles:
        pf.print_message("", "info")
        pf.print_message("Creating Blast database with complete subject schema...", "info")
        blastdb_path: str = os.path.join(blast_folder, 'blastdb')
        ff.create_directory(blastdb_path)
        blast_db_files: str = os.path.join(blastdb_path, 'genbank_protein_db')
        makeblastdb_exec: str = lf.get_tool_path('makeblastdb')
        bf.make_blast_db(makeblastdb_exec, master_file_path, blast_db_files, 'prot')

        # Run BLAST rep query vs alleles subject
        pf.print_message("Running BLASTs (rep vs allele) between schemas...", "info")
        best_bsr_values: Dict[str, Dict[str, float]] = run_blasts_match_schemas(query_translations_rep_paths,
                                                                                blast_db_files,
                                                                                blast_folder,
                                                                                self_score_dict,
                                                                                max_id_length,
                                                                                get_blastp_exec,
                                                                                bsr,
                                                                                cpu)
        pf.print_message("", None)

        # Remove matched id from the full subject schema file
        locus_removal = 0
        subject_base_list = set()
        for query, match in best_bsr_values.items():
            for subject, score in match.items():
                subject_base = re.sub(r'_\*?\d+$', '', subject)
                subject_base_list.add(subject_base)
        
        for subject in subject_base_list:
            subject_translations_paths.pop(subject, None)
            locus_removal += 1

        # Write the best BLAST matches to a file
        pf.print_message("Writting results to the output file...", "info")
        best_blast_matches_file = write_best_blast_matches_to_file(best_bsr_values, query_translations_rep_paths, subject_translations_paths, output_directory, rep_vs_alleles, 'rep_vs_alleles')

        pf.print_message(f"From the rep vs alleles Blast {len(best_bsr_values)} matches were found.", "info")
        pf.print_message(f"{locus_removal} loci were removed.", "info")
        pf.print_message(f"{len(subject_translations_paths)} subject loci had no matches.", "info")

    # Clean up temporary files
    if not no_cleanup:
        pf.print_message("Cleaning up temporary files...", "info")
        # Remove temporary files
        ff.cleanup(output_directory, [best_blast_matches_file, logf.get_log_file_path(gb.LOGGER)])

    return best_blast_matches_file