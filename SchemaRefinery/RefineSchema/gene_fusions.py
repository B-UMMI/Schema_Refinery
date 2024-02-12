import os
import sys
import pickle
import concurrent.futures
from itertools import repeat

try:
    from utils import (file_functions as ff,
                       sequence_functions as sf,
                       clustering_functions as cf,
                       blast_functions as bf,
                       alignments_functions as af,
                       kmers_functions as kf,
                       list_functions as lf)
except ModuleNotFoundError:
    from SchemaRefinery.utils import (file_functions as ff,
                                      sequence_functions as sf,
                                      clustering_functions as cf,
                                      blast_functions as bf,
                                      alignments_functions as af,
                                      kmers_functions as kf,
                                      list_functions as lf)

def fetch_not_included_cds(file_path_cds):
    """
    Compares the hashes list with the hashes obtained from cds.

    Parameters
    ----------
    hash_list : list
        List containing all of the sequences hashes present in the input files.

    Returns
    -------
    not_included_cds : dict
        Returns dict with key as fasta header and value as fasta sequence.
    """
    
    not_included_cds = {}
    i = 1
    # Read FASTA files
    for rec in sf.read_fasta_file_iterator(file_path_cds):
        print(f"Processed {i} CDS")
        i += 1
        not_included_cds[rec.id] = rec.seq

    return not_included_cds


def alignment_dict_to_file(blast_results_dict, file_path):
    """
    Writes alignments strings to file.

    Parameters
    ----------
    alignment_string_dict : dict
        dict containing alignments string to write to file
    file_path : str
        File path to create to write the file

    Returns
    -------
    No return
    """
    # Write first column into TSV file
    with open(file_path, 'w') as report_file:
        report_file.writelines(["Query\t",
                                "Subject\t",
                                "Query length\t",
                                "Subject length\t",
                                "Query start\t",
                                "Query end\t",
                                "Subject start\t",
                                "Subject end\t",
                                "Length\t",
                                "Score\t",
                                "Number of Gaps\t",
                                "Pident - Percentage of identical matches\t",
                                "prot_BSR\t",
                                "prot_seq_Kmer_sim\t",
                                "prot_seq_Kmer_cov\t",
                                "cluster_frequency_in_genomes_query_cds\t",
                                "cluster_frequency_in_genomes_subject_cds\t",
                                "palign_global\t",
                                "palign_local\t",
                                "class\n"])
        # Write all of the strings into TSV file
        for results in blast_results_dict.values():
            for result in results.values():
                for r in result.values():
                    report_file.writelines('\t'.join([str(r) for r in r.values()]) + '\n')


def separate_blastn_results_into_classes(representative_blast_results, path):
    """
    Separates one BLASTn dict into various classes and adds them into one dict

    Parameters
    ----------
    representative_blast_results : dict
        dict containing BLAST results
    path : str
        Dir path to write the files for each class of results

    Returns
    -------
    cluster_classes : dict
        dict containing cluster separated by classes
    """

    def add_to_class_dict(class_name):
        """
        Function that adds entries to the cluster_classes dict based on class.
        
        Parameters
        ----------
        class_name : str
            Class name, which is a key in cluster_classes to where to add entries.
            
        Returns
        -------
        Modifies the dict created in the parent function.
        """
        
        cluster_classes[class_name][query][id_subject].update({id_: blastn_entry})
        
    def add_class_to_dict(class_name):
        """
        Adds class as last item in representative_blast_results dict
        
        Parameters
        ----------
        class_name : str
            Class name, which is a key in cluster_classes to where to add entries.
            
        Returns
        -------
        Modifies the representative_blast_results dict in the parent function.
        """
        
        representative_blast_results[query][id_subject][id_].update({'class': class_name})
        
        
    # Create various dicts and split the results into different classes
    cluster_classes = {}
    results_outcome = {}
    classes = ['1a',
               '1b',
               '1c',
               '2',
               '3',
               '4']
    # Add all of the classes keys and their query keys into the dict
    for class_ in classes:
        cluster_classes[class_] = {}
        # create dict for each rep inside cluster_classes
        for query in representative_blast_results.keys():
            cluster_classes[class_][query] = {}
            for subjects in representative_blast_results.keys():
                cluster_classes[class_][query][subjects] = {}

    results_outcome['Join'] = []
    results_outcome['Retain'] = set()    
    # Process results into classes
    for query, rep_blast_result in representative_blast_results.items():
        for id_subject, matches in rep_blast_result.items():
            for id_, blastn_entry in matches.items():
                
                if blastn_entry['global_palign'] >= 0.8:
                    # Based on BSR
                    if blastn_entry['bsr'] >= 0.6:
                        add_class_to_dict('1a')
                        add_to_class_dict('1a')
                        results_outcome['Join'].append([query, id_subject])
                    # If BSR <0.6 verify if query cluster is the most prevalent
                    elif blastn_entry['frequency_in_genomes_query_cds'] >= blastn_entry['frequency_in_genomes_subject_cds'] * 10:
                        add_class_to_dict('1b')
                        add_to_class_dict('1b')
                        results_outcome['Retain'].add(query)
                    # Add two as separate
                    else:
                        add_class_to_dict('1c')
                        add_to_class_dict('1c')
                        results_outcome['Retain'].add(query)
                        results_outcome['Retain'].add(id_subject)
                # Palign < 0.8        
                else:
                    # verify if query cluster is the most prevalent
                    if blastn_entry['frequency_in_genomes_query_cds'] >= blastn_entry['frequency_in_genomes_subject_cds'] * 10:
                        add_class_to_dict('3')
                        add_to_class_dict('3')
                        results_outcome['Retain'].add(query)
                    else:
                        add_class_to_dict('2')
                        add_to_class_dict('2')
    
    # Join the various CDS groups into single group based on ids matches    
    results_outcome['Join'] = [join for join in cf.cluster_by_ids(results_outcome['Join'])]
    
    # Remove blastn results that are present in another classes thus removing
    # empty query entries
    for class_id, entries in list(cluster_classes.items()):
        for query_id, results in list(entries.items()):
            for subjects in list(results):
                if len(cluster_classes[class_id][query_id][subjects]) == 0:
                    del cluster_classes[class_id][query_id][subjects]
                    
            if len(cluster_classes[class_id][query_id]) == 0:
                del cluster_classes[class_id][query_id]
                
    # Write all of the classes into TSV files
    for k, v in cluster_classes.items():
        report_file_path = os.path.join(path, f"blastn_group_{k}.tsv")
        # Write individual class to file
        alignment_dict_to_file(v, report_file_path)

    return cluster_classes, results_outcome

def decode_CDS_sequences_ids(path_to_file):
    """
    Function to read a dict contained in pickle file and decode its values based
    on polyline.
    
    Parameters
    ----------
    path_to_file : str
        Path to the pickle file.
        
    Returns
    -------
    decoded_dict : dict
        Contains hashed CDS as keys, and number id of the genome which
        that CDS is present
    
    """
    
    with open(path_to_file, "rb") as infile:
        hash_table = pickle.load(infile)
    
    decoded_dict = {}
    for key, value in hash_table.items():
        decoded_dict[key] = lf.polyline_decoding(value)
        
    return decoded_dict
    
def main(schema, output_directory, allelecall_directory, clustering_sim,
         clustering_cov, alignment_ratio_threshold_gene_fusions,
         pident_threshold_gene_fusions, genome_presence, cpu):

    temp_folder = os.path.join(allelecall_directory,"temp")
    file_path_cds = os.path.join(allelecall_directory, "unclassified_sequences.fasta")
    
    if not os.path.exists(temp_folder) or not os.path.exists(file_path_cds):
        sys.exit(f"Error: {temp_folder} must exist, make sure that AlleleCall "
                 "was run using --no-cleanup and --output-unclassified flag.")

    # Verify if the dataset is small, if it is, keep minimum genomes in which
    # specific CDS cluster is present to 5 if not to 1% of the dataset size.
    if not genome_presence:
        count_genomes_path = os.path.join(temp_folder, "1_cds_prediction")
        number_of_genomes = len(os.listdir(count_genomes_path))
        if number_of_genomes <= 20:
            genome_presence = 5
        else:
            genome_presence = round(number_of_genomes * 0.01)
            
    # Put all constants in one dict in order to decrease number of variables
    # used around.
    constants_threshold = [alignment_ratio_threshold_gene_fusions, 
                           pident_threshold_gene_fusions,
                           genome_presence]
        
    print("Identifying CDS present in the schema...")
    cds_present = os.path.join(temp_folder,"3_cds_preprocess/cds_deduplication/distinct_cds_merged.hashtable")
    
    decoded_sequences_ids = decode_CDS_sequences_ids(cds_present)

    print("Identifying CDS not present in the schema...")
        
    not_included_cds = fetch_not_included_cds(file_path_cds)

    print(f"Identified {len(not_included_cds)} valid CDS not present in the schema")
    
    ff.create_directory(output_directory)
    cds_output = os.path.join(output_directory, "1_CDS_processing")
    ff.create_directory(cds_output)
    cds_not_present_file_path = os.path.join(cds_output, "CDS_not_found.fasta")
    
    frequency_cds = {}
    with open(cds_not_present_file_path, 'w+') as cds_not_found:
        for id_, sequence in not_included_cds.items():
            cds_not_found.writelines(">"+id_+"\n")
            cds_not_found.writelines(str(sequence)+"\n")
            
            hashed_seq = sf.seq_to_hash(str(sequence)) 
            if hashed_seq in decoded_sequences_ids:
                frequency_cds[id_] = len(decoded_sequences_ids[hashed_seq]) - 1
            else:
                frequency_cds[id_] = 0
                

    print("Translate unclassified CDS...")
    cds_translation_dict = {}
    cds_not_present_translation_file_path = os.path.join(cds_output, "CDS_not_found_translation.fasta")
    protein_hashes = {}
    i = 1
    total = len(not_included_cds)
    with open(cds_not_present_translation_file_path, 'w+') as translation:
        for id_s, sequence in not_included_cds.items():
            print(f"Translated {i}/{total} CDS")
            i += 1
            protein_translation = str(sf.translate_dna(str(sequence),
                                                       11,
                                                       0,
                                                       True)[0][0])
            
            prot_hash = sf.seq_to_hash(protein_translation)
            if prot_hash not in protein_hashes:
                protein_hashes[prot_hash] = [id_s]
                cds_translation_dict[id_s] = protein_translation
                translation.writelines(id_s+"\n")
                translation.writelines(protein_translation+"\n")
            else:
                protein_hashes[prot_hash].append(id_s)

    print("Extracting minimizers for the translated sequences and clustering...")

    reps_groups = {}
    clusters = {}
    reps_sequences = {}

    # sort by size of proteins
    cds_translation_dict = {k: v for k, v in sorted(cds_translation_dict.items(),
                                                    key=lambda x: len(x[1]),
                                                    reverse=True)}

    [clusters, reps_sequences, 
     reps_groups, prot_len_dict] = cf.minimizer_clustering(cds_translation_dict,
                                                           5, 5, True, 1, 
                                                           clusters,
                                                           reps_sequences, 
                                                           reps_groups,
                                                           1, clustering_sim, 
                                                           clustering_cov,
                                                           True)
    print("Filtering clusters...")
    # Get frequency of cluster
    frequency_cds_cluster = {rep: sum([frequency_cds[entry[0]] for entry in value]) for rep, value in clusters.items()}
    # Filter cluster by the total sum of CDS that are present in the genomes, based on input value
    clusters = {rep: cluster_member for rep, cluster_member in clusters.items() if frequency_cds_cluster[rep] >= constants_threshold[2]}

    print("Retrieving kmers similiarity and coverage between representatives...")
    reps_kmers_sim = {}
    reps_translation_dict = {rep_id: rep_seq for rep_id, rep_seq in cds_translation_dict.items()
                             if rep_id in clusters}

    reps_translation_dict = {k: v for k, v in sorted(reps_translation_dict.items(),
                                                     key=lambda x: len(x[1]),
                                                     reverse=True)}
    # recalculate the sim and cov between reps, get all of the values, so threshold
    # is set to 0.
    for cluster_id in reps_translation_dict:
        kmers_rep = set(kf.determine_minimizers(reps_translation_dict[cluster_id],
                                                5, 5, 1, True, True))
        
        reps_kmers_sim[cluster_id] = cf.select_representatives(kmers_rep,
                                                               reps_groups,
                                                               0,
                                                               0,
                                                               prot_len_dict,
                                                               cluster_id, 5)

        reps_kmers_sim[cluster_id] = {match_values[0]: match_values[1:]
                                      for match_values in reps_kmers_sim[cluster_id]}

    print("Running BLASTn between cluster representatives...")
    blast_output = os.path.join(output_directory, "2_BLAST_processing")
    ff.create_directory(blast_output)
    
    blastn_output = os.path.join(blast_output, "BLASTn_processing")
    ff.create_directory(blastn_output)
    # write fasta file with all of the representatives sequences
    representatives_blastn_folder = os.path.join(blastn_output,
                                                "cluster_representatives_fastas")
    ff.create_directory(representatives_blastn_folder)

    representatives_all_fasta_file = os.path.join(representatives_blastn_folder,
                                                  "all_cluster_representatives.fasta")

    rep_paths_nuc = {}
    # Write FASTA files
    with open(representatives_all_fasta_file, 'w') as all_fasta:
        for cluster_rep_id in clusters:

            all_fasta.writelines(">"+cluster_rep_id+"\n")
            all_fasta.writelines(str(not_included_cds[cluster_rep_id])+"\n")

            rep_fasta_file = os.path.join(representatives_blastn_folder,
                                          f"cluster_rep_{cluster_rep_id}.fasta")
            rep_paths_nuc[cluster_rep_id] = rep_fasta_file
            with open(rep_fasta_file, 'w') as rep_fasta:
                rep_fasta.writelines(">"+cluster_rep_id+"\n")
                rep_fasta.writelines(str(not_included_cds[cluster_rep_id])+"\n")

    blastn_results_folder = os.path.join(blastn_output, "blast_results")
    ff.create_directory(blastn_results_folder)

    total_reps = len(rep_paths_nuc)
    representative_blast_results = {}
    representative_blast_results_coords = {}
    i = 1
    # Run BLASTn for all representatives
    with concurrent.futures.ProcessPoolExecutor(max_workers=cpu) as executor:
        for res in executor.map(bf.run_all_representative_blasts_multiprocessing,
                                clusters.keys(),
                                repeat('blastn'),
                                repeat(blastn_results_folder),
                                repeat(rep_paths_nuc),
                                repeat(representatives_all_fasta_file)):

            filtered_alignments_dict, _, alignment_coords = af.get_alignments_dict_from_blast_results(
                res[1], constants_threshold[1], True)
            # Save the BLASTn results
            representative_blast_results.update(filtered_alignments_dict)
            representative_blast_results_coords.update(alignment_coords)

            print(
                f"Running BLASTn for cluster representatives: {res[0]} - {i}/{total_reps}")
            i += 1

    print("Running BLASTp based on BLASTn results matches...")
    # Obtain the list for what BLASTp runs to do, no need to do all vs all as previously
    blastp_runs_to_do = {query: lf.flatten_list([[query],[subject[1]['subject']
                                            for subject in subjects.values()]]) 
                         for query, subjects in representative_blast_results.items()}
    
    # Create all of the folders
    blastp_results = os.path.join(blast_output,
                                  "BLASTp_processing")
    ff.create_directory(blastp_results)
    
    blastn_results_matches_translations = os.path.join(blastp_results,
                                                       "blastn_results_matches_translations")
    ff.create_directory(blastn_results_matches_translations)

    representatives_blastp_folder = os.path.join(blastn_results_matches_translations,
                                                "cluster_rep_translation")
    ff.create_directory(representatives_blastp_folder)
    

    
    blastp_results_folder = os.path.join(blastp_results,
                                         "BLASTp_results")
    ff.create_directory(blastp_results_folder)
    
    rep_paths_prot = {}
    rep_matches_prot = {}    
    # Write the protein FASTA files
    for query_id, subjects_ids in blastp_runs_to_do.items():
        rep_translation_file = os.path.join(representatives_blastp_folder,
                                            f"cluster_rep_translation_{query_id}.fasta")
        rep_paths_prot[query_id] = rep_translation_file
        with open(rep_translation_file, 'w') as trans_fasta:
            trans_fasta.writelines(">"+query_id+"\n")
            trans_fasta.writelines(str(reps_translation_dict[query_id])+"\n")
            
        rep_matches_translation_file = os.path.join(blastn_results_matches_translations,
                                                    f"cluster_matches_translation_{query_id}.fasta")
        
        rep_matches_prot[query_id] = rep_matches_translation_file
        with open(rep_matches_translation_file, 'w') as trans_fasta:            
            for subject_id in subjects_ids:
                trans_fasta.writelines(">"+subject_id+"\n")
                trans_fasta.writelines(str(reps_translation_dict[subject_id])+"\n")


    # Calculate BSR based on BLASTp.
    total_blasts = len(blastp_runs_to_do)
    bsr_values = {}
    # Create query entries
    for query in blastp_runs_to_do.keys():
        bsr_values[query] = {}
        
    i = 1
    with concurrent.futures.ProcessPoolExecutor(max_workers=cpu) as executor:
        for res in executor.map(bf.run_all_representative_blasts_multiprocessing,
                                blastp_runs_to_do.keys(), 
                                repeat('blastp'),
                                repeat(blastp_results_folder),
                                repeat(rep_paths_prot),
                                rep_matches_prot.values()):
            
            filtered_alignments_dict, self_score, _ = af.get_alignments_dict_from_blast_results(res[1], 0, False)
            # Since BLAST may find several local aligments choose the largest one to calculate BSR
            for query, subjects_dict in filtered_alignments_dict.items():
                for subject_id, results in subjects_dict.items():
                    largest_alignment = 0
                    for entry_id, result in results.items():
                        if result['query_length'] > largest_alignment:
                            largest_alignment = result['query_length']
                            bsr_values[query].update({subject_id: bf.compute_bsr(result['score'], self_score)})

            print(f"Running BLASTp for cluster representatives matches: {res[0]} - {i}/{total_blasts}")
            i += 1

    # Add kmer cov, kmer sim and frequency of the cds in the genomes
    for query, subjects_dict in list(representative_blast_results.items()):
        for subject, blastn_results in subjects_dict.items():
            # Some results may not have kmers matches or BSR values so put them
            # as 0
            if subject in reps_kmers_sim[query]:
                sim, cov = reps_kmers_sim[query][subject]
            else:
                sim = 0
                cov = 0
            if subject in bsr_values[query]:
                bsr = bsr_values[query][subject]
                # For some reason some isolates have bsr slighty higher than 1
                # related to blast database and sequences used.
                if bsr > 1.0:
                    bsr = round(bsr)
            else:
                bsr = 0
                
            # Calculate total alignment for all of the fragments of BLASTn
            # if there more than one BLASTn matches
            # For query and subject
            if len(blastn_results) > 1:
                total_length = {}
                gaps = representative_blast_results_coords[query][subject].pop('gaps')
                for ref, intervals in representative_blast_results_coords[query][subject].items():
                    # Sort by start position
                    sorted_intervals = [interval for interval in sorted(intervals,
                                                                        key=lambda x: x[0])]
                    # Merge alignments and calculate the total length of BLASTn alignments.
                    length = sum([(interval[1] - interval[0] + 1) for interval in af.merge_intervals(sorted_intervals)])
                    total_length[ref] = length
                # Calculate global palign
                global_palign = min([(total_length['query'] + gaps) / blastn_results[1]['query_length'],
                                     (total_length['subject']) / blastn_results[1]['subject_length']])
            # If there is only on results
            else:
                global_palign = min([(blastn_results[1]['query_end'] - blastn_results[1]['query_start'] + 1 + result['gaps']) / blastn_results[1]['query_length'],
                                    (blastn_results[1]['subject_end'] - blastn_results[1]['subject_start'] + 1) / blastn_results[1]['subject_length']])                
            # Create update dict with the values to add
            update_dict = {'bsr' : bsr,
                           'kmers_sim': sim,
                           'kmers_cov': cov,
                           'frequency_in_genomes_query_cds' : frequency_cds_cluster[query],
                           'frequency_in_genomes_subject_cds' : frequency_cds_cluster[subject],
                           'global_palign': global_palign}
            
            for entry_id, result in list(blastn_results.items()):
                # Calculate Palign
                local_palign = min([(result['query_end'] - result['query_start'] + 1 + result['gaps']) / result['query_length'],
                                    (result['subject_end'] - result['subject_start'] + 1) / result['subject_length']])
                
                if local_palign >= 0:
                    # update the update dict
                    update_dict.update({'local_palign' : local_palign})
                    # Add everything to the dict
                    representative_blast_results[query][subject][entry_id].update(update_dict)
                # Remove palign that is negative, meaning tha reverse blast match was made
                else:
                    del representative_blast_results[query][subject][entry_id]

    print("Filtering BLASTn results into subclusters...")
    blastn_processed_results = os.path.join(blast_output, "Processed_Blastn")
    ff.create_directory(blastn_processed_results)
    report_file_path = os.path.join(blastn_processed_results, "blastn_all_matches.tsv")
    
    # Separate results into different classes
    cluster_classes, results_outcome = separate_blastn_results_into_classes(representative_blast_results,
                                                                            blastn_processed_results)
    # Write all of the BLASTn results to a file
    alignment_dict_to_file(representative_blast_results, report_file_path)