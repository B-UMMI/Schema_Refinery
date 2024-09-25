import os
import concurrent.futures
import itertools
import statistics
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
    

def identify_paralagous_loci(schema_directory, output_directory, cpu_cores, blast_score_ratio,
                             translation_table, size_threshold, run_mode):
    # Identify all of the fastas in the schema directory
    fasta_files_dict = {loci.split('.')[0]: os.path.join(schema_directory, loci) for loci in os.listdir(schema_directory) if os.path.isfile(os.path.join(schema_directory, loci))}
    # Identify all of the fastas short in the schema directory
    short_folder = os.path.join(schema_directory, 'short')
    fasta_files_short_dict = {loci.split('.')[0]: os.path.join(short_folder, loci) for loci in os.listdir(short_folder) if os.path.isfile(os.path.join(short_folder, loci))}
    
    blast_folder = os.path.join(output_directory, 'Blast')
    ff.create_directory(blast_folder)
    translation_folder = os.path.join(output_directory, 'Translation')
    ff.create_directory(translation_folder)
    len_short_folder = len(fasta_files_dict)
    master_file_path = os.path.join(blast_folder, 'master_file.fasta')
    query_paths_dict = {}
    protein_size_mode_dict = {}
    i = 1
    for loci in fasta_files_dict:
        subject_fasta = fasta_files_dict[loci] if run_mode == 'alleles_vs_alleles' or run_mode == 'reps_vs_alleles' else fasta_files_short_dict[loci]
        query_fasta = fasta_files_dict[loci] if run_mode == 'alleles_vs_alleles' else fasta_files_short_dict[loci]
        query_fasta_translation = os.path.join(translation_folder, f"{loci}-translation.fasta")
        query_paths_dict[loci] = query_fasta_translation

        print(f"\rTranslated loci FASTA: {i}/{len_short_folder}", end='', flush=True)
        i += 1
        fasta_dict = sf.fetch_fasta_dict(query_fasta, False)
        with open(query_fasta_translation, 'w') as query_file:
            loci_allele_size = []
            for allele_id, sequence in fasta_dict.items():
                loci_allele_size.append(len(sequence))
                protseq = sf.translate_sequence(str(sequence), translation_table)
                query_file.write(f">{allele_id}\n{str(protseq)}\n")

        # Calculate the mode of the lengths of the sequences
        if loci_allele_size:
            protein_size_mode_dict[loci] = statistics.mode(loci_allele_size)

        fasta_dict = sf.fetch_fasta_dict(subject_fasta, False)
        write_type = 'w' if not os.path.exists(master_file_path) else 'a'
        with open(master_file_path, write_type) as master_file:
            for allele_id, sequence in fasta_dict.items():
                protseq = sf.translate_sequence(str(sequence), translation_table)
                master_file.write(f">{allele_id}\n{str(protseq)}\n")


    # For better prints
    max_id_length = len(max(query_paths_dict.keys(), key=len))
    # Blast executable
    blast_exec = lf.get_tool_path('blastp')
    # Self-score folder
    blastp_results_ss_folder = os.path.join(blast_folder, 'blastp_results_ss')
    ff.create_directory(blastp_results_ss_folder)

    self_score_dict = {}
    i = 1
    # Calculate self-score
    print("\nCalculating self-score for each loci:")
    with concurrent.futures.ProcessPoolExecutor(max_workers=cpu_cores) as executor:
        for res in executor.map(bf.run_self_score_multiprocessing,
                                query_paths_dict.keys(),
                                itertools.repeat(blast_exec),
                                query_paths_dict.values(),
                                itertools.repeat(blastp_results_ss_folder)):
            
            _, self_score, _, _ = af.get_alignments_dict_from_blast_results(res[1], 0, False, True, True, True, True)
    
            # Save self-score
            self_score_dict[res[0]] = self_score
                            
            print(f"\rRunning BLASTp to calculate self-score for {res[0]: <{max_id_length}}", end='', flush=True)
            i += 1    
    
    # Print newline
    print('\n')  

    # Get makeblastdb executable
    makeblastdb_exec = lf.get_tool_path('makeblastdb')
    # Path to the blast database folder
    blast_db = os.path.join(blast_folder, 'Blast_db_prot')
    ff.create_directory(blast_db)
    # Path to the blast database files
    blast_db_prot = os.path.join(blast_db, 'Blast_db_protein')
    bf.make_blast_db(makeblastdb_exec, master_file_path, blast_db_prot, 'prot')
    # Blast output folder
    blast_output_folder = os.path.join(blast_folder, 'Blast_output')
    ff.create_directory(blast_output_folder)

    bsr_values = {}
    best_bsr_values = {}
    total_blasts = len(query_paths_dict)
    i = 1
    print(f"\nRunning BLASTp for each loci {'alleles' if run_mode == 'alleles_vs_alleles' else 'representatives'}"
          f" vs {'representatives' if run_mode == 'reps_vs_reps' else 'alleles'}:")
    with concurrent.futures.ProcessPoolExecutor(max_workers=cpu_cores) as executor:
        for res in executor.map(bf.run_blastdb_multiprocessing, 
                                itertools.repeat(blast_exec),
                                itertools.repeat(blast_db_prot),
                                query_paths_dict.values(),
                                query_paths_dict.keys(),
                                itertools.repeat(blast_output_folder)):
            
            filtered_alignments_dict, _, _, _ = af.get_alignments_dict_from_blast_results(res[1], 0, True, False, True, True, False)

            # Since BLAST may find several local aligments choose the largest one to calculate BSR.
            for query, subjects_dict in filtered_alignments_dict.items():
                query_loci_id = query.split('_')[0]
                best_bsr_values.setdefault(query_loci_id, {})
                bsr_values.setdefault(query, {})
                for subject_id, results in subjects_dict.items():
                    #Highest score (First one)
                    subject_score = next(iter(results.values()))['score']
                    computed_score = bf.compute_bsr(subject_score, self_score_dict[res[0]])
                    bsr_values[query].update({subject_id: computed_score})
                    # Get the best BSR value between loci
                    # If the BSR is better than the current best, update it
                    subject_loci_id = subject_id.split('_')[0]
                    if not best_bsr_values[query_loci_id].get(subject_loci_id) or best_bsr_values[query_loci_id][subject_loci_id] >= computed_score:
                        best_bsr_values[query_loci_id][subject_loci_id] = computed_score
                    
                    
            print(f"\rRunning BLASTp for cluster representatives matches: {res[0]} - {i}/{total_blasts: <{max_id_length}}", end='', flush=True)
            i += 1
            
    # Print newline
    print('\n')
    
    filtered_best_bsr_values = {
    query_loci_id: {
        subject_loci_id: computed_score
        for subject_loci_id, computed_score in subject_dict.items()
        if computed_score >= blast_score_ratio
    }
    for query_loci_id, subject_dict in best_bsr_values.items()
    }
    
    paralagous_loci_report = os.path.join(output_directory, 'paralagous_loci_report.tsv')
    paralagous_list = []
    paralagous_list_mode_check = []
    with open(paralagous_loci_report, 'w') as report_file:
        report_file.write("Query_loci_id\tSubject_loci_id\tBSR\tMode_check\n")
        for query_loci_id, subject_dict in filtered_best_bsr_values.items():
            for subject_loci_id, computed_score in subject_dict.items():
                paralagous_list.append((query_loci_id, subject_loci_id))
                mode_check = stats.modes_within_value(protein_size_mode_dict[query_loci_id], protein_size_mode_dict[subject_loci_id], size_threshold)
                report_file.write(f"{query_loci_id}\t{subject_loci_id}\t{computed_score}\t{mode_check}\n")
                
                if mode_check:
                    paralagous_list_mode_check.append((query_loci_id, subject_loci_id))
                
    
    paralagous_list = cf.cluster_by_ids(paralagous_list)
    paralagous_loci_report_cluster_by_id = os.path.join(output_directory, 'paralagous_loci_report_cluster_by_id.tsv')
    with open(paralagous_loci_report_cluster_by_id, 'a') as report_file:
        for cluster in paralagous_list:
            report_file.write(f"{cluster}\n")
            
    paralagous_list_mode_check = cf.cluster_by_ids(paralagous_list_mode_check)
    paralagous_loci_report_mode = os.path.join(output_directory, 'paralagous_loci_report_mode.tsv')
    with open(paralagous_loci_report_mode, 'a') as report_file:
        for cluster in paralagous_list_mode_check:
            report_file.write(f"{cluster}\n")