import os
import statistics
import sys
import pathlib
import shutil
from Bio import SeqIO
from typing import Dict, List, Tuple

try:
    from SchemaAnnotation import (consolidate as cs)
    from utils import (sequence_functions as sf,
                        blast_functions as bf,
                        linux_functions as lf,
                        file_functions as ff,
                        alignments_functions as af,
                        statistics as stats,
                        clustering_functions as cf,
                        print_functions as pf,
                        logger_functions as logf,
                        constants as ct,
                        globals as gb)
except ModuleNotFoundError:
    from SchemaRefinery.SchemaAnnotation import (consolidate as cs)
    from SchemaRefinery.utils import (sequence_functions as sf,
                                    blast_functions as bf,
                                    linux_functions as lf,
                                    file_functions as ff,
                                    alignments_functions as af,
                                    statistics as stats,
                                    clustering_functions as cf,
                                    print_functions as pf,
                                    logger_functions as logf,
                                    constants as ct,
                                    globals as gb)



def identify_paralogous_loci(schema_directory: str, 
                             output_directory: str,
                             annotation_paths: List[str], 
                             cpu: int, 
                             bsr: float,
                             translation_table: int, 
                             size_threshold: float, 
                             processing_mode: str,
                             no_cleanup: bool,) -> None:
    """
    Identify paralogous loci by performing BLAST searches and analyzing sequence similarities.

    Parameters
    ----------
    schema_directory : str
        Path to the directory containing schema FASTA files.
    output_directory : str
        Path to the directory where output files will be saved.
    annotation_paths : List[str]
        Paths to the files with annotations to be added to the final recommendations file.
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

    pf.print_message("Identify Schema FASTAs...", "info")

    # Identify all of the fastas in the schema directory
    fasta_files_dict: Dict[str, str] = {
        loci.split('.')[0]: os.path.join(schema_directory, loci)
        for loci in os.listdir(schema_directory)
        if os.path.isfile(os.path.join(schema_directory, loci)) and loci.endswith('.fasta')
    }
    # Identify all of the fastas short in the schema directory
    short_folder: str = os.path.join(schema_directory, 'short')
    fasta_files_short_dict: Dict[str, str] = {
        loci.split('.')[0].rsplit('_', 1)[0]: os.path.join(short_folder, loci)
        for loci in os.listdir(short_folder)
        if os.path.isfile(os.path.join(short_folder, loci)) and loci.endswith('.fasta')
    }

    output_d= os.path.abspath(output_directory)
    blast_folder: str = os.path.join(output_d, 'Blast')
    ff.create_directory(blast_folder)
    translation_folder: str = os.path.join(blast_folder, 'Translation')
    ff.create_directory(translation_folder)
    # Total loci
    len_short_folder: int = len(fasta_files_dict)
    master_file_path: str = os.path.join(blast_folder, 'master_file.fasta')
    query_paths_dict: Dict[str, str] = {}
    all_loci_allele_size_stats = {}
    i: int = 1
    pf.print_message('')
    pf.print_message('Translating fasta files', 'info')
    # Translate the sequences that are to be used in the BLASTp search
    for loci in fasta_files_dict:
        # Get the subject and query FASTA files
        subject_fasta: str = fasta_files_dict[loci] if processing_mode.split('_')[-1] == 'alleles' else fasta_files_short_dict[loci]
        query_fasta: str = fasta_files_short_dict[loci] if processing_mode.split('_')[0] == 'reps' else fasta_files_dict[loci]
        fetch_size_fasta: str = fasta_files_dict[loci]
        # Translation variables
        query_fasta_translation: str = os.path.join(translation_folder, f"{loci}_translation.fasta")
        query_paths_dict[loci] = query_fasta_translation

        pf.print_message(f"Translated loci FASTA: {i}/{len_short_folder}", "info", end='\r', flush=True)
        i += 1
        # Get the sizes for each loci (all alleles)
        loci_allele_size: List[int] = []
        fasta_dict: Dict[str, str] = sf.fetch_fasta_dict(fetch_size_fasta, False)
        for allele_id, sequence in fasta_dict.items():
            loci_allele_size.append(len(sequence))

        # Calculate the mode of the lengths of the sequences
        if loci_allele_size:
            all_loci_allele_size_stats[loci] = (int(min(loci_allele_size)/3),
                                                int(max(loci_allele_size)/3),
                                                int(statistics.mode(loci_allele_size)/3),
                                                statistics.mean(loci_allele_size)/3)
    
        # Get the fasta sequences for the query
        fasta_dict = sf.fetch_fasta_dict(query_fasta, False)
        # Write the sequences to the query file
        with open(query_fasta_translation, 'w') as query_file:
            for allele_id, sequence in fasta_dict.items():
                protseq: str = sf.translate_sequence(str(sequence), translation_table)
                query_file.write(f">{allele_id}\n{str(protseq)}\n")

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

    pf.print_message('')
    pf.print_message('Calculating self-scores', 'info')
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
    # Run BLASTp in parallel
    pf.print_message('')
    pf.print_message('Running Blastp', 'info')
    blastp_results_files = bf.run_blastp_operations(cpu,
                                                    blast_exec,
                                                    blast_db_prot,
                                                    query_paths_dict,
                                                    blast_output_folder,
                                                    total_blasts,
                                                    max_id_length)
    
    for blast_result_file in blastp_results_files:
        # Get the filtered alignments
        filtered_alignments_dict, _, _, _ = af.get_alignments_dict_from_blast_results(blast_result_file,
                                                                                        0,
                                                                                        True,
                                                                                        False,
                                                                                        True,
                                                                                        True,
                                                                                        False)

        # Since BLAST may find several local alignments choose the largest one to calculate BSR.
        for query, subjects_dict in filtered_alignments_dict.items():
            # Only use the alleles that have a self_score
            # Not all alleles pass the threshoolds when calculating the self-score
            if query in self_score_dict.keys():
                query_loci_id: str = query.rsplit('_', 1)[0]
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
                    # Since there may be several paralogous we save the various matches
                    subject_loci_id: str = subject_id.rsplit('_', 1)[0]
                    if not best_bsr_values[query_loci_id].get(subject_loci_id):
                        best_bsr_values[query_loci_id][subject_loci_id] = computed_score
                    elif computed_score > best_bsr_values[query_loci_id][subject_loci_id]:
                        best_bsr_values[query_loci_id][subject_loci_id] = computed_score

    pf.print_message(f"", None)
    
    paralogous_loci_report: str = os.path.join(output_d, 'paralogous_loci_report.tsv')
    paralogous_list: List[Tuple[str, str]] = []
    paralogous_list_check: List[Tuple[str, str]] = []
    # Write the report file with all of the paralogous loci results
    pf.print_message('Writting report file with all of the paralogous loci results', 'info')
    with open(paralogous_loci_report, 'w') as report_file:
        report_file.write("Query_loci_id\t"
                          "Subject_loci_id\t"
                          "BSR\t"
                          "if_loci_intersect\t"
                          "if_close_distance\t"
                          "Loci_min_allele_size\t"
                          "Loci_max_allele_size\t"
                          "Loci_mode_allele_size\t"
                          "Loci_mean_allele_size\t"
                          "\n")
        # Check if the loci are paralogous
        for query_loci_id, subject_dict in best_bsr_values.items():
            for subject_loci_id, computed_score in subject_dict.items():
                paralogous_list.append((query_loci_id, subject_loci_id))

                if_loci_intersect = stats.if_loci_intersect(all_loci_allele_size_stats[query_loci_id][:2],
                                                        all_loci_allele_size_stats[subject_loci_id][:2])
                
                if not if_loci_intersect:
                    if_close_distance = stats.calculate_loci_distance(all_loci_allele_size_stats[query_loci_id][:3],
                                                                    all_loci_allele_size_stats[subject_loci_id][:3],
                                                                    size_threshold)
                else:
                    if_close_distance = True

                report_file.write(f"{query_loci_id}\t"
                                  f"{subject_loci_id}\t"
                                  f"{computed_score}\t"
                                  f"{if_loci_intersect}\t"
                                  f"{if_close_distance}\t"
                                  f"{all_loci_allele_size_stats[query_loci_id][0]}|{all_loci_allele_size_stats[subject_loci_id][0]}\t"
                                  f"{all_loci_allele_size_stats[query_loci_id][1]}|{all_loci_allele_size_stats[subject_loci_id][1]}\t"
                                  f"{all_loci_allele_size_stats[query_loci_id][2]}|{all_loci_allele_size_stats[subject_loci_id][2]}\t"
                                  f"{all_loci_allele_size_stats[query_loci_id][3]}|{all_loci_allele_size_stats[subject_loci_id][3]}\t"
                                  "\n")
                
                if if_loci_intersect or if_close_distance:
                    paralogous_list_check.append((query_loci_id, subject_loci_id))


    # open the sample file used 
    paralogous = open(paralogous_loci_report)
    if len(paralogous.readlines()) < 2:
        pf.print_message('No paralogs were found', 'info')
        # Clean up temporary files
        if not no_cleanup:
            pf.print_message("Cleaning up temporary files...", "info")
            # Remove temporary files
            ff.cleanup(output_d, [logf.get_log_file_path(gb.LOGGER)])
        sys.exit(0)

    # Write recomendation file
    pf.print_message('Writting recommendations TSV file', 'info')
    header: str = "Locus\tAction\n"
    Joined: List[str] = []
    par_loci = 0
    paralogous_list_check = cf.cluster_by_ids(paralogous_list_check)
    paralogous_loci_report_mode: str = os.path.join(output_d, 'paralogous_loci_final_recommendations.tsv')
    with open(paralogous_loci_report_mode, 'w') as report_file:
        report_file.write(header)
        for cluster in paralogous_list_check:
            for loci in cluster:
                # The action is always 'Join'
                report_file.write(f"{loci}\tJoin\n")
                Joined.append(loci)
                par_loci += 1
            # Each cluster is separated by a '#' row
            report_file.write("#\t\n")
        for loci, loci_path in query_paths_dict.items():
            if loci not in Joined:
               report_file.write(f"{loci}\tAdd\n")

    # Cluster the paralogous loci by id and write the results to a file
    pf.print_message('Writting Clusters by IDs TSV file', 'info')
    header: str = "Joined_loci_id\tClustered_loci_ids\n"
    paralogous_list = cf.cluster_by_ids(paralogous_list)
    paralogous_loci_report_cluster_by_id: str = os.path.join(output_d, 'paralogous_loci_report_cluster_by_id.tsv')
    with open(paralogous_loci_report_cluster_by_id, 'w') as report_file:
        report_file.write(header)
        for cluster in paralogous_list:
            report_file.write(f"Joined_{cluster[0]}\t{','.join(cluster)}\n#\n")

    # Final statistics
    pf.print_message('')
    pf.print_message(f'{par_loci} loci were found to be paralogous and joined into {len(paralogous_list_check)} groups.', 'info')
    pf.print_message('')

    # Annotate the recomendations files using the consolidate option from the annotation module
    consolidated_annotations = os.path.join(output_d, "paralogous_annotations.tsv") 
    if annotation_paths:
        # Append the recommendation files to the begginig of the list of files to annotate
        files: List[str]
        files = [paralogous_loci_report_mode] + annotation_paths
        pf.print_message("Consolidating annoations...", "info")
        # Cleanup is put as False by deafult
        consolidated_annotations: str = cs.consolidate_annotations(files,
                                    False,
                                    consolidated_annotations)
         # Final statistics
        annotations_count = 0
        hypoteticals = 0
        with open(consolidated_annotations, "r") as f:
            for line in f:
                parts = line.strip().split("\t")
                if len(parts) > 3:
                    if any(field != "NA" for field in parts[1:]):
                        annotations_count += 1
                if any("hypothetical protein" in field.lower() for field in parts):
                    hypoteticals += 1

        
        prf.print_message(f'A total of {annotations_count-1} loci were annotated.', 'info')
        prf.print_message(f'From these {hypoteticals} loci where annotated as "hypothetical proteins".', 'info')

        pf.print_message('Annotation consolidation successfully completed.', 'info')
        pf.print_message('')


    # Clean up temporary files
    if not no_cleanup:
        pf.print_message("Cleaning up temporary files...", "info")
        # Remove temporary files
        ff.cleanup(output_d, [paralogous_loci_report,
                                    paralogous_loci_report_mode,
                                    consolidated_annotations if annotation_paths else None,
                                    logf.get_log_file_path(gb.LOGGER)])