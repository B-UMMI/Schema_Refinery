import os
import concurrent.futures
import itertools
import statistics
from typing import Dict, List, Tuple

try:
    from utils import (
                        sequence_functions as sf,
                        blast_functions as bf,
                        linux_functions as lf,
                        file_functions as ff,
                        alignments_functions as af,
                        statistics as stats,
                        clustering_functions as cf,
    )
except ModuleNotFoundError:
    from SchemaRefinery.utils import (
                                    sequence_functions as sf,
                                    blast_functions as bf,
                                    linux_functions as lf,
                                    file_functions as ff,
                                    alignments_functions as af,
                                    statistics as stats,
                                    clustering_functions as cf,
    )

def identify_paralogous_loci(schema_directory: str, 
                             output_directory: str, 
                             cpu: int, 
                             bsr: float,
                             translation_table: int, 
                             size_threshold: float, 
                             processing_mode: str) -> None:
    """
    Identify paralogous loci by performing BLAST searches and analyzing sequence similarities.

    Parameters
    ----------
    schema_directory : str
        Path to the directory containing schema FASTA files.
    output_directory : str
        Path to the directory where output files will be saved.
    cpu : int
        Number of CPU cores to use for parallel processing.
    bsr : float
        Threshold for the BLAST score ratio to consider loci as paralogous.
    translation_table : int
        Translation table number to use for translating nucleotide sequences to protein sequences.
    size_threshold : float
        Threshold for the size difference to consider loci as paralogous.
    processing_mode: str
        Processing mode to determine which sequences to use for BLAST.

    Returns
    -------
    None

    Notes
    -----
    - The function creates several output files in the specified output directory.
    - It performs BLAST searches to identify paralogous loci and writes the results to output files.
    - The function handles different processing modes to determine which sequences to use for BLAST.

    Examples
    --------
    >>> schema_directory = '/path/to/schema'
    >>> output_directory = '/path/to/output'
    >>> cpu = 4
    >>> bsr = 0.8
    >>> translation_table = 11
    >>> size_threshold = 0.2
    >>> identify_paralogous_loci(schema_directory, output_directory, cpu, bsr, translation_table, size_threshold, processing_mode)
    """

    # Identify all of the fastas in the schema directory
    fasta_files_dict: Dict[str, str] = {
        loci.split('.')[0]: os.path.join(schema_directory, loci)
        for loci in os.listdir(schema_directory)
        if os.path.isfile(os.path.join(schema_directory, loci)) and loci.endswith('.fasta')
    }
    # Identify all of the fastas short in the schema directory
    short_folder: str = os.path.join(schema_directory, 'short')
    fasta_files_short_dict: Dict[str, str] = {
        loci.split('.')[0].split('_')[0]: os.path.join(short_folder, loci)
        for loci in os.listdir(short_folder)
        if os.path.isfile(os.path.join(short_folder, loci)) and loci.endswith('.fasta')
    }
    blast_folder: str = os.path.join(output_directory, 'Blast')
    ff.create_directory(blast_folder)
    translation_folder: str = os.path.join(output_directory, 'Translation')
    ff.create_directory(translation_folder)
    len_short_folder: int = len(fasta_files_dict)
    master_file_path: str = os.path.join(blast_folder, 'master_file.fasta')
    query_paths_dict: Dict[str, str] = {}
    protein_size_mode_dict: Dict[str, int] = {}
    i: int = 1
    # Translate the sequences that are to be used in the BLASTp search
    for loci in fasta_files_dict:
        # Get the subject and query FASTA files
        subject_fasta: str = fasta_files_dict[loci] if processing_mode.split('_')[-1] == 'alleles' else fasta_files_short_dict[loci]
        query_fasta: str = fasta_files_short_dict[loci] if processing_mode.split('_')[0] == 'rep' else fasta_files_dict[loci]
        query_fasta_translation: str = os.path.join(translation_folder, f"{loci}-translation.fasta")
        query_paths_dict[loci] = query_fasta_translation

        print(f"\rTranslated loci FASTA: {i}/{len_short_folder}", end='', flush=True)
        i += 1
        # Get the fasta sequences for the query
        fasta_dict: Dict[str, str] = sf.fetch_fasta_dict(query_fasta, False)
        # Write the sequences to the query file
        with open(query_fasta_translation, 'w') as query_file:
            loci_allele_size: List[int] = []
            for allele_id, sequence in fasta_dict.items():
                loci_allele_size.append(len(sequence))
                protseq: str = sf.translate_sequence(str(sequence), translation_table)
                query_file.write(f">{allele_id}\n{str(protseq)}\n")

        # Calculate the mode of the lengths of the sequences
        if loci_allele_size:
            protein_size_mode_dict[loci] = statistics.mode(loci_allele_size)
        # Get the fasta sequences for the subject
        fasta_dict = sf.fetch_fasta_dict(subject_fasta, False)
        # Write the sequences to the master file
        write_type: str = 'w' if not os.path.exists(master_file_path) else 'a'
        with open(master_file_path, write_type) as master_file:
            for allele_id, sequence in fasta_dict.items():
                protseq = sf.translate_sequence(str(sequence), translation_table)
                master_file.write(f">{allele_id}\n{str(protseq)}\n")

    # For better prints
    max_id_length: int = len(max(query_paths_dict.keys(), key=len))
    # Blast executable
    blast_exec: str = lf.get_tool_path('blastp')
    # Calculate self-score
    self_score_dict: Dict[str, float] = bf.calculate_self_score(query_paths_dict,
                                                                blast_exec,
                                                                blast_folder,
                                                                max_id_length,
                                                                cpu)   

    # Get makeblastdb executable
    makeblastdb_exec: str = lf.get_tool_path('makeblastdb')
    # Path to the blast database folder
    blast_db: str = os.path.join(blast_folder, 'Blast_db_prot')
    ff.create_directory(blast_db)
    # Path to the blast database files
    blast_db_prot: str = os.path.join(blast_db, 'Blast_db_protein')
    bf.make_blast_db(makeblastdb_exec, master_file_path, blast_db_prot, 'prot')
    # Blast output folder
    blast_output_folder: str = os.path.join(blast_folder, 'Blast_output')
    ff.create_directory(blast_output_folder)

    bsr_values: Dict[str, Dict[str, float]] = {}
    best_bsr_values: Dict[str, Dict[str, float]] = {}
    total_blasts: int = len(query_paths_dict)
    i = 1
    print(f"\nRunning BLASTp...")
    with concurrent.futures.ProcessPoolExecutor(max_workers=cpu) as executor:
        for res in executor.map(bf.run_blastdb_multiprocessing, 
                                itertools.repeat(blast_exec),
                                itertools.repeat(blast_db_prot),
                                query_paths_dict.values(),
                                query_paths_dict.keys(),
                                itertools.repeat(blast_output_folder)):
            
            filtered_alignments_dict, _, _, _ = af.get_alignments_dict_from_blast_results(res[1], 0, True, False, True, True, False)

            # Since BLAST may find several local alignments choose the largest one to calculate BSR.
            for query, subjects_dict in filtered_alignments_dict.items():
                query_loci_id: str = query.split('_')[0]
                best_bsr_values.setdefault(query_loci_id, {})
                bsr_values.setdefault(query, {})
                for subject_id, results in subjects_dict.items():
                    # Highest score (First one)
                    subject_score: float = next(iter(results.values()))['score']
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
                    # Get the best BSR value between loci
                    # If the BSR is better than the current best, update it
                    subject_loci_id: str = subject_id.split('_')[0]
                    if not best_bsr_values[query_loci_id].get(subject_loci_id):
                        best_bsr_values[query_loci_id][subject_loci_id] = computed_score
                    elif computed_score > best_bsr_values[query_loci_id][subject_loci_id]:
                        best_bsr_values[query_loci_id][subject_loci_id] = computed_score
                    
            print(f"\rRunning BLASTp for cluster representatives matches: {res[0]} - {i}/{total_blasts: <{max_id_length}}", end='', flush=True)
            i += 1
            
    # Print newline
    print('\n')
    
    paralogous_loci_report: str = os.path.join(output_directory, 'paralogous_loci_report.tsv')
    paralogous_list: List[Tuple[str, str]] = []
    paralogous_list_mode_check: List[Tuple[str, str]] = []
    # Write the report file with all of the paralogous loci results
    with open(paralogous_loci_report, 'w') as report_file:
        report_file.write("Query_loci_id\tSubject_loci_id\tBSR\tMode_check\n")
        for query_loci_id, subject_dict in best_bsr_values.items():
            for subject_loci_id, computed_score in subject_dict.items():
                paralogous_list.append((query_loci_id, subject_loci_id))
                mode_check: bool = stats.modes_within_value(protein_size_mode_dict[query_loci_id], protein_size_mode_dict[subject_loci_id], size_threshold)
                report_file.write(f"{query_loci_id}\t{subject_loci_id}\t{computed_score}\t{mode_check}\n")
                
                if mode_check:
                    paralogous_list_mode_check.append((query_loci_id, subject_loci_id))
    # Cluster the paralogous loci by id and write the results to a file
    paralogous_list = cf.cluster_by_ids(paralogous_list)
    paralogous_loci_report_cluster_by_id: str = os.path.join(output_directory, 'paralogous_loci_report_cluster_by_id.tsv')
    with open(paralogous_loci_report_cluster_by_id, 'a') as report_file:
        for cluster in paralogous_list:
            report_file.write(f"{','.join(cluster)}\n")
    # Cluster the paralogous loci by id that passed the mode check and write the results to a file
    paralogous_list_mode_check = cf.cluster_by_ids(paralogous_list_mode_check)
    paralogous_loci_report_mode: str = os.path.join(output_directory, 'paralogous_loci_report_mode.tsv')
    with open(paralogous_loci_report_mode, 'a') as report_file:
        for cluster in paralogous_list_mode_check:
            report_file.write(f"Joined_{cluster[0]}\t{','.join(cluster)}\n#\n")
