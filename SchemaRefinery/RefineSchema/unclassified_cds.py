import os
import concurrent.futures
from itertools import repeat

try:
    from utils import (file_functions as ff,
                       sequence_functions as sf,
                       clustering_functions as cf,
                       blast_functions as bf,
                       alignments_functions as af,
                       kmers_functions as kf,
                       iterable_functions as itf,
                       graphical_functions as gf)
except ModuleNotFoundError:
    from SchemaRefinery.utils import (file_functions as ff,
                                      sequence_functions as sf,
                                      clustering_functions as cf,
                                      blast_functions as bf,
                                      alignments_functions as af,
                                      kmers_functions as kf,
                                      iterable_functions as itf,
                                      graphical_functions as gf)

def alignment_dict_to_file(blast_results_dict, file_path, write_type):
    """
    Writes alignments strings to file.

    Parameters
    ----------
    alignment_string_dict : dict
        dict containing alignments string to write to file
    file_path : str
        File path to create to write the file
    write_type : str
        If to create new file and write or to append to existing file

    Returns
    -------
    No return, writes or appends a file at the file_path
    """
    
    header = ["Query\t",
              "Subject\t",
              "Query_length\t",
              "Subject_length\t",
              "Query_start\t",
              "Query_end\t",
              "Subject_start\t",
              "Subject_end\t",
              "Length\t",
              "Score\t",
              "Number_of_gaps\t",
              "Pident\t",
              "Prot_BSR\t",
              "Prot_seq_Kmer_sim\t",
              "Prot_seq_Kmer_cov\t",
              "Cluster_frequency_in_genomes_query_cds\t",
              "Cluster_frequency_in_genomes_subject_cds\t",
              "Global_palign_all\t",
              "Global_palign_pident_min\t",
              "Global_palign_pident_max\t",
              "Palign_local\t",
              "Class\n"]
    
    # Write or append to the file
    with open(file_path, write_type) as report_file:
        # Write the header only if the file is being created
        if write_type == 'w':
            report_file.write("".join(header))
        # Write all the alignment data
        for results in blast_results_dict.values():
            for result in results.values():
                for r in result.values():
                    report_file.write('\t'.join(map(str, r.values())) + '\n')

def add_items_to_results(representative_blast_results, reps_kmers_sim, bsr_values,
                         representative_blast_results_coords_all,
                         representative_blast_results_coords_pident,
                         frequency_cds_cluster):
    """
    Function to add to BLAST results additional information, it adds:
    bsr: value between the two CDS.
    kmers_sim: kmer similarities.
    kmer_cov: kmer coverage.
    frequency_in_genomes_query_cds: how many times that query appears in the schema genomes.
    frequency_in_genomes_subject_cds: how many subject appers in the schema genomes.
    global_palign_all: minimum of how much query or subject covers each other.
    global_palign_pident_min: minimum of how much query or subject covers each other, takes into
    account only entries with specific pident value.
    global_palign_pident_max: maximum of how much query or subject covers each other, takes into
    account only entries with specific pident value.
                                  
    Parameters
    ----------                                  
    representative_blast_results : 
        Dict that contains representatibes BLAST results.
    reps_kmers_sim : dict
        Dict that contains values for kmer similarities between CDS.
    bsr_values : dict
        Dict that contains BSR values between CDS.
    representative_blast_results_coords_all : dict
        Dict that contain the coords for all of the entries.
    representative_blast_results_coords_pident : dict
        Dict that contain the coords for all of the entries above certain pident value.
    frequency_cds_cluster : dict
        Dict that contains sum of frequency of that representatives cluster in the
        genomes of the schema.
    
    Returns
    -------
    No returns, modifies the representative_blast_results dict inside the main
    function.
    """
    # Add kmer cov, kmer sim and frequency of the cds in the genomes
    for query, subjects_dict in list(representative_blast_results.items()):
        for subject, blastn_results in list(subjects_dict.items()):
            # Some results may not have kmers matches or BSR values so put them
            # as 0
            if subject in reps_kmers_sim[query]:
                sim, cov = reps_kmers_sim[query][subject]
            else:
                sim = 0
                cov = 0
            
            # Get BSR value, if not in query return 0, meaning that even though
            # there was BLASTn match no BLASTp match was made.
            bsr = bsr_values[query].get(subject, 0)
            # For some reason some isolates have bsr slighty higher than 1
            # related to blast database and sequences used.
            if bsr > 1.0:
                bsr = float(round(bsr))
            # Calculate total alignment for all of the fragments of BLASTn
            # if there more than one BLASTn alignments
            # For query and subject
            for entry_id, result in list(blastn_results.items()):
                total_length = {}
                # Sum all the intervals
                for ref, intervals in representative_blast_results_coords_all[query][subject].items():
                    sorted_intervals = sorted(intervals, key=lambda x: x[0])
                    length = sum(interval[1] - interval[0] + 1 for interval in af.merge_intervals(sorted_intervals))
                    total_length[ref] = length
                # Calculate global palign
                global_palign_all = min(total_length['query'] / result['query_length'],
                                        total_length['subject'] / result['subject_length'])
                # Sum the intervals with desired pident threshold
                for ref, intervals in representative_blast_results_coords_pident[query][subject].items():
                    if intervals:
                        sorted_intervals = sorted(intervals, key=lambda x: x[0])
                        length = sum(interval[1] - interval[0] + 1 for interval in af.merge_intervals(sorted_intervals))
                        total_length[ref] = length
                    else:
                        total_length[ref] = 0
                
                # Calculated min and max palign
                global_palign_pident_min = min(total_length['query'] / result['query_length'],
                                                total_length['subject'] / result['subject_length'])

                global_palign_pident_max = max(total_length['query'] / result['query_length'],
                                                total_length['subject'] / result['subject_length'])
                # Calculate local palign
                local_palign = min((result['query_end'] - result['query_start'] + 1) / result['query_length'],
                                    (result['subject_end'] - result['subject_start'] + 1) / result['subject_length'])
                # If the alignment is more than 0
                if local_palign >= 0:
                    update_dict = {
                        'bsr': bsr,
                        'kmers_sim': sim,
                        'kmers_cov': cov,
                        'frequency_in_genomes_query_cds': frequency_cds_cluster[query],
                        'frequency_in_genomes_subject_cds': frequency_cds_cluster[subject],
                        'global_palign_all': global_palign_all,
                        'global_palign_pident_min': global_palign_pident_min,
                        'global_palign_pident_max': global_palign_pident_max,
                        'local_palign': local_palign
                    }
                    representative_blast_results[query][subject][entry_id].update(update_dict)
                # If not then the inverse alignment was made, so we remove it
                else:
                    del representative_blast_results[query][subject][entry_id]
                    
            if not representative_blast_results[query][subject]:
                del representative_blast_results[query][subject]
        if not representative_blast_results[query]:
            del representative_blast_results[query]

def separate_blastn_results_into_classes(representative_blast_results, constants):
    """
    Separates one BLASTn dict into various classes and adds them as entry to the
    representative_blast_results dict

    Parameters
    ----------
    representative_blast_results : dict
        Dict that contains representatibes BLAST results with all of the additional
        info.
    constants : list
        Contains the constants to be used in this function
        
    Returns
    -------
    classes_outcome : list
        List of list that contains class IDS used in the next function
    drop_list : set
        Contains the CDS IDs to be removed from further processing to appearing
        fewer time in genomes than their match.
    """
        
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
    
    # Create all of the classes order based on priority when choosing which CDS
    # to keep and remove e.g if two CDS are joined by one while another result makes
    # them separete we keep them joined or when one result X and Y chooses class
    # 2a which means that X is 10x more frequent in the genomes while another
    # results Y and Z adds them as separate in this case we only keep X and Z
    classes_outcome = ['1a',
                       '1b',
                       '2a',
                       '1c',
                       '2b',
                       '3a',
                       '3b']
    # Set of IDs to drop
    drop_list = set()
    # Process results into classes
    for query, rep_blast_result in representative_blast_results.items():
        for id_subject, matches in rep_blast_result.items():
            for id_, blastn_entry in matches.items():
                query_subject_freq = blastn_entry['frequency_in_genomes_query_cds']/blastn_entry['frequency_in_genomes_subject_cds']
                subject_query_freq = blastn_entry['frequency_in_genomes_subject_cds']/blastn_entry['frequency_in_genomes_query_cds']
                
                if blastn_entry['global_palign_all'] >= 0.8:
                    # Based on BSR
                    if blastn_entry['bsr'] >= 0.6:
                        add_class_to_dict('1a')
                    # If BSR <0.6 verify if between CDS there ir more than 10x
                    # difference in presence in the schema
                    elif min([query_subject_freq,subject_query_freq]) >= 0.1:
                        add_class_to_dict('1b')
                        # Remove the subject
                        if blastn_entry['frequency_in_genomes_query_cds'] > blastn_entry['frequency_in_genomes_query_cds']:
                            drop_list.add(id_subject)
                        # Remove the query
                        else:
                            drop_list.add(query)
                    # Add two as separate
                    else:
                        add_class_to_dict('1c')
                # Palign < 0.8        
                else:
                    if blastn_entry['pident'] >= constants[1]:
                        # Verify if between CDS there ir more than 10x
                        # difference in presence in the schema
                        if min([query_subject_freq,subject_query_freq]) >= 0.1:
                            add_class_to_dict('2a')
                            # Remove the subject
                            if blastn_entry['frequency_in_genomes_query_cds'] > blastn_entry['frequency_in_genomes_query_cds']:
                                drop_list.add(id_subject)
                            # Remove the query
                            else:
                                drop_list.add(query)
                        # If in similiar proportion in genomes then add two as separate
                        else:
                            add_class_to_dict('2b')
                            
                    else:
                        # If one CDS is contained inside another
                        if blastn_entry['global_palign_pident_max'] >= 0.8:
                            add_class_to_dict('3a')
                        # Everything else not classified
                        else:
                            add_class_to_dict('3b')
    
    return classes_outcome, drop_list

def process_classes(representative_blast_results, classes_outcome, drop_list):
    """
    Process the classified representative_blast_results to identify the CDS
    that are to be kept as potential new loci while also adding the different
    relationships between these CDS for further processing.
    
    Parameters
    ----------
    representative_blast_results : dict
        Dict that contains representatibes BLAST results with all of the additional
        info.
    classes_outcome : list
        All of the existing classes.
    drop_list
        Contains the CDS IDs to be removed from further processing to appearing
        fewer time in genomes than their match.
        
    Returns
    -------
    cds_to_keep : dict     
        Dict of the CDS to keep by each classification.
    relationships : dict
        Dict that contains relationships between various CDS and clusters.
    """
    # Create variables.
    # Variable to add the CDS what will be kept by class.
    cds_to_keep = {}
    # Variable for other ralationships between CDS e.g if two CDS X and Y are 
    # classified as 1a (join) and there are other classfication with them
    # for example Z and X classified as 2b (keep both) then we add to this dict
    # this case for end-user to know.
    other_relationships = {}
    # List that contains the various CDS with 1a classification to join.
    cluster_to_join = []
    
    for class_ in classes_outcome:
        cds_to_keep[class_] = set()
        other_relationships[class_] = []
        # Get all entries with the desired class.
        class_dict = {query : {subject: {id_: entry for id_, entry in entries.items() if entry['class'] == class_}
                                   for subject, entries in subjects.items()}
                          for query, subjects in representative_blast_results.items()}
        # Remove empty dicts entries
        itf.remove_empty_dicts_recursive(class_dict)
        
        # Process the CDS to find what CDS to retain while also adding the
        # relationships between different CDS.
        for query, rep_blast_result in class_dict.items():
            for id_subject, matches in rep_blast_result.items():
                
                # Find cases that were already processed or to be dropped.
                processed_cases = itf.flatten_list([[c for c in cds] for cds in cds_to_keep.values()])
                processed_cases.append(drop_list)
                # Process all of the cases that have 1a classification.
                # even if they may be in drop_list
                if class_ == '1a':
                    cds_to_keep[class_].update([query, id_subject])
                    cluster_to_join.append([query, id_subject])
                    # Since query and subject CDS are kept inside joined clusters
                    # then remove them from drop_list
                    drop_list = [x for x in drop_list if x not in [query, id_subject]]
                    
                # All of the other classifications
                else:
                    # Get those cases that query and subject were not processed.
                    if query not in processed_cases and id_subject not in processed_cases:
                        cds_to_keep[class_].update([query, id_subject])
                    # If query was not processed.
                    elif query not in processed_cases:
                        cds_to_keep[class_].add(query)
                    # If subject was not processed.
                    elif id_subject not in processed_cases:
                        cds_to_keep[class_].add(id_subject)
                    # If to add relationship between different CDS and clusters.
                    # Keep only 3a as they are the one relevant for the output
                    if (query not in drop_list or id_subject not in drop_list) and class_ == '3a':
                        other_relationships[class_].append([query, id_subject])

    # Create the joined cluster by joining by IDs.
    cds_to_keep['1a'] = {i+1: join for i, join in enumerate(cf.cluster_by_ids(cluster_to_join))}
            
    # And remove duplicates
    for class_, results in other_relationships.items():
        other_relationships[class_] = itf.get_unique_sublists(other_relationships[class_])
        
    # Initialize relationships dictionary for each class.
    relationships = {class_: {} for class_ in other_relationships}
    
    # Process all relationships identified (Identify the Joined cluster ID).
    for class_, results in other_relationships.items():
        # Classification 3b is to be ignored because the is no need for them
        # it means that query and subject are too different.
        if class_ == '3b':
            continue
        # Iterate over results
        for result in results:
             # Identify the cluster ID, if it is joined id or CDS ID.
            cluster_id = itf.identify_string_in_dict(result[1], cds_to_keep['1a']) or result[1]
            # Update relationships dictionary.
            relationships[class_].setdefault(cluster_id, []).append(result)

    # Get unique relationships.
    for class_, relationship in list(relationships.items()):
        for id_, r in list(relationship.items()):
            relationships[class_][id_] = itf.get_unique_sublists(r)
            
    return cds_to_keep, relationships

def write_processed_results_to_file(cds_to_keep, relationships, representative_blast_results,
                                    classes_outcome, output_path):
    """
    Write the results from processed_classes into various files.
    
    Parameters
    ----------
    cds_to_keep : dict
        Dict of the CDS to keep by each classification.
    relationships : dict
        Dict that contains relationships between various clusters.
    representative_blast_results : dict
        Dict that contains representatibes BLAST results with all of the additional
        info.
    classes_outcome : list
        List of list that contains class IDS used in the next function.
    output_path : str
        Path were to write files.
        
    Returns
    -------
    No returns, writes files in output path.
    """
    # Create directory
    blast_by_cluster_output = os.path.join(output_path, "blast_results_by_cluster")
    ff.create_directory(blast_by_cluster_output)
    
    # Get all of the BLAST entries for that cluster.
    for class_, cds in cds_to_keep.items():
        for i, cluster in enumerate(cds):
            if class_ == '1a':
                id_ = i + 1
                cluster_type = 'joined_cluster'
                cluster = cds[id_]
            else:
                id_ = cluster
                cluster = [cluster]
                cluster_type = 'retained'
            write_dict = {query : {subject: {id_: entry for id_, entry in entries.items()}
                                   for subject, entries in subjects.items()}
                          for query, subjects in representative_blast_results.items()
                          if query in cluster}
        
            report_file_path = os.path.join(blast_by_cluster_output, f"blast_{cluster_type}_{id_}.tsv")
            alignment_dict_to_file(write_dict, report_file_path, 'w')
    
    # Create directory.
    joined_cluster_relationships_output = os.path.join(blast_by_cluster_output, "1_blast_results_by_cluster_relationships")
    ff.create_directory(joined_cluster_relationships_output)
    report_relationships_output = os.path.join(blast_by_cluster_output, "2_relationships_to_joined_clusters")
    ff.create_directory(report_relationships_output)
    
    # Write blast results by class and relationships.
    for class_, relationship in relationships.items():
        # Skip empty entries.
        if not relationship:
            continue
        # Process and write relationships.
        for cluster_id, r_ids in relationship.items():
            # Get all IDs into a list.
            query_ids = [query_id[0] for query_id in r_ids]
            subject_ids = [subject_id[1] for subject_id in r_ids]
            
            # Get transformed values into the list (replaces CDS name for joined cluster ID).
            transformed_query_ids = [itf.identify_string_in_dict(id_, cds_to_keep['1a']) or id_ for id_ in query_ids]
            
            # Get entries based on BLAST.
            write_dict = {query : {subject: {id_: entry for id_, entry in entries.items()
                                             if entry['class'] == class_}
                                   for subject, entries in subjects.items() if subject in subject_ids}
                          for query, subjects in representative_blast_results.items() if query in query_ids}
                
            report_file_path = os.path.join(joined_cluster_relationships_output, f"blast_relationships_to_{cluster_id}.tsv")
            # What write type to use (if file already exists).
            write_type = 'a' if os.path.exists(report_file_path) else 'w'
            # Write BLAST results to file.
            alignment_dict_to_file(write_dict, report_file_path, write_type)
            # Path to the relationships report.
            relationships_report_file_path = os.path.join(report_relationships_output, f"relationships_to_cluster_{cluster_id}_report.txt")
            
            # Write all of the report files.
            with open(relationships_report_file_path, write_type) as relationships_report_file:
                # Class that has CDS that are partially contained or cantains other CDS.
                if class_ == '3a':

                    relationships_report_file.writelines("The following CDS may partially contain the following "
                                                         "CDS from this cluster:\n")
                else:
                    relationships_report_file.writelines("The following CDS matched with BLASTn to the following"
                                                         f" elements of this cluster and have the classification '{class_}'"
                                                         " to this elements however they are probably different loci:\n")
                # Write the cluster ID.
                relationships_report_file.writelines(f"{cluster_id}:\n")
                seen = ""
                for i, query_id in enumerate(transformed_query_ids):
                    # To increase redability and reduce redundancy we simplify the output.
                    if query_id == seen:
                        white_spaces = itf.create_whitespace_string(f"{'CDS' if type(query_id) == str else 'Cluster'} {query_id} entry: {query_ids[i]} against ")
                        relationships_report_file.write("\t" + white_spaces + f"{subject_ids[i]}\n")
                    else:
                        relationships_report_file.write(f"\t{'CDS' if type(query_id) == str else 'Cluster'} {query_id} entry: {query_ids[i]} against {subject_ids[i]}\n")
                    
                    if type(query_id) == str:
                        seen = query_id
                    else:
                        seen = query_ids[i]
    # Write all of the ids inside Joined cluster.
    # Create directory .
    cluster_members_output = os.path.join(output_path, "joined_cluster_members")
    ff.create_directory(cluster_members_output)
    # Write files.
    for cluster_id, cluster in cds_to_keep['1a'].items():
        cluster_output_path = os.path.join(cluster_members_output, f"Joined_cluster_{cluster_id}.txt")
        with open(cluster_output_path, 'w') as output:
            for c in cluster:
                output.writelines(c + '\n')

    # Create directory.
    joined_cluster_relationships_output = os.path.join(output_path, "blast_results_by_class")
    ff.create_directory(joined_cluster_relationships_output)
    # Write classes to file.
    for class_ in classes_outcome:
        # Fetch all entries with the desired class.
        write_dict = {query : {subject: {id_: entry for id_, entry in entries.items() if entry['class'] == class_}
                               for subject, entries in subjects.items()}
                      for query, subjects in representative_blast_results.items()}

        report_file_path = os.path.join(joined_cluster_relationships_output, f"blastn_group_{class_}.tsv")
        # Write individual class to file.
        alignment_dict_to_file(write_dict, report_file_path, 'w')
            
def wrap_up_blast_results(cds_to_keep, not_included_cds, clusters,
                          output_path, constants, cpu):
    """
    This function wraps up the results for this module by writing FASTAs files
    for the possible new loci to include into the schema and creates graphs for
    each results group. It also translates schema short FASTAs into proteins
    and the possible new loci, to calculate the BSR values to see if those
    possible new loci are already present in the schema.
    
    Parameters
    ----------
    cds_to_keep : dict
        Dict of the CDS to keep by each classification.
    not_included_cds : dict
        Dict that contains all of the DNA sequences for all of the CDS.
    clusters : dict
        Dict that contains the cluster representatives as keys and similar CDS
        as values.
    output_path : str
        Path to were write the FASTA files.
    cpu : int
        Number of cores to use in BLAST multiprocessing.
        
    Returns
    -------
    outcomes_trans_reps : dict
        Contains the translations of the representatives.
    """
    # Create directories.
    cds_outcome_results_graphs = os.path.join(output_path, "Blast_results_outcomes_graphs")
    ff.create_directory(cds_outcome_results_graphs)
    
    fasta_folder = os.path.join(output_path, "results_outcomes_fastas")
    ff.create_directory(fasta_folder)
    
    cds_outcome_results_fastas_folder = os.path.join(fasta_folder, "results_outcomes_fastas")
    ff.create_directory(cds_outcome_results_fastas_folder)
    
    cds_outcome_results_reps_fastas_folder = os.path.join(fasta_folder, "results_outcomes_reps_fastas")
    ff.create_directory(cds_outcome_results_reps_fastas_folder)
    
    # Write FASTA files for each CDS group to join or retain.
    print("Writting FASTA file for possible new loci...")
    outcome_paths = {}
    outcome_paths_reps = {}
    for class_, cds_list in cds_to_keep.items():
        i = 1
        for cds in cds_list:
            if class_ == '1a':
                class_name_cds = f"joined_{i}"
            elif class_ == '3a':
                class_name_cds = f"for_reference_{class_}_{cds}"  
            elif class_ == '3b':
                class_name_cds = f"thrash_{class_}_{cds}"
            else:
                class_name_cds = f"retained_{class_}_{cds}"
            
            # Create the files paths
            cds_group_fasta_file = os.path.join(cds_outcome_results_fastas_folder, class_name_cds + '.fasta')
            outcome_paths[class_name_cds] = cds_group_fasta_file
            
            cds_group_reps_file = os.path.join(cds_outcome_results_reps_fastas_folder, class_name_cds + '.fasta')
            outcome_paths_reps[class_name_cds] = cds_group_reps_file
            # to decrease the number of code lines just iterate over string id as a list
            if type(cds) == str:
                cds = [cds]
            else:
                cds = cds_to_keep[class_][cds]
            # Add to counter
            i += 1
            # Write all of the CDS in the group + all of the sequences in the
            # cluster of those representatives
            with open(cds_group_fasta_file, 'w') as fasta_file:
                for rep_id in cds:
                    cds_ids = [cds_id for cds_id in clusters[rep_id]]
                    for cds_id in cds_ids:
                        fasta_file.writelines(f">{cds_id}\n")
                        fasta_file.writelines(str(not_included_cds[cds_id])+"\n")
            # Write the representatives CDS
            with open(cds_group_reps_file, 'w') as fasta_file:
                for rep_id in cds:
                    fasta_file.writelines(f">{rep_id}\n")
                    fasta_file.writelines(str(not_included_cds[rep_id])+"\n")
                    
    # Create directories.
    cds_outcome_trans = os.path.join(fasta_folder, "cds_outcome_translation_deduplicated")
    ff.create_directory(cds_outcome_trans)
    # Translate possible new loci.
    outcomes_trans = {}
    for key, o_path in outcome_paths.items():
        trans_path = os.path.join(cds_outcome_trans, key + ".fasta")
        outcomes_trans[key] = trans_path
        fasta_dict = sf.fetch_fasta_dict(o_path, False)
        trans_dict, _, _ = sf.translate_seq_deduplicate(fasta_dict, trans_path,
                                                        None,
                                                        constants[5], 
                                                        False)

    # Translate all clusters into their respective cluster file.
    # Create directories
    cds_outcome_trans_rep = os.path.join(fasta_folder, "cds_outcome_translation_reps_deduplicated")
    ff.create_directory(cds_outcome_trans_rep)
    # Translate possible new loci representatives.
    outcomes_trans_reps = {}
    for key, o_path in outcome_paths_reps.items():
        trans_path = os.path.join(cds_outcome_trans, key + ".fasta")
        outcomes_trans_reps[key] = trans_path
        fasta_dict = sf.fetch_fasta_dict(o_path, False)
        trans_dict, _, _ = sf.translate_seq_deduplicate(fasta_dict, 
                                                        trans_path, 
                                                        None,
                                                        constants[5], 
                                                        False)

    # Create graphs for all results and for each class.
    classes_tsv_path = os.path.join(output_path, 'blast_results_by_class')
    for tsv_file_path in os.listdir(classes_tsv_path):
        abs_path = os.path.join(classes_tsv_path, tsv_file_path)
        file_name = os.path.basename(tsv_file_path)
        # Create directories.
        graphs_path = os.path.join(cds_outcome_results_graphs, f"{file_name.replace('tsv','')}")
        ff.create_directory(graphs_path)
        # Render histograms.
        gf.render_histogram(abs_path,
                            graphs_path,
                            ['Query_length', 'Subject_length'],
                            ['Length', 'Count'])
        # Render line charts.
        gf.render_line_chart(abs_path,
                             graphs_path,
                             ['Pident', 'Prot_BSR', 'Prot_seq_Kmer_sim',
                              'Prot_seq_Kmer_cov'],
                             ['Entries', 'Values'], False)
    return outcomes_trans_reps

def process_schema(schema, outcomes_translations_reps, output_path, self_score_dict,
                   constants, cpu):
    """
    This function processes data related to the schema seed, importing, translating
    and BLASTing against the unclassified CDS clusters representatives to validate
    them.
    
    Parameters
    ----------
    schema : str
        Path to the schema seed folder.
    outcomes_translations_reps : dict
        Dict with the paths to the translations of the unclassified CDS clusters.
    output_path : str
        Path were to write the results of this function.
    self_score : dict
        Self-score for BSR calculation of the unclassified CDS reps to consider.
    cpu : int
        Number of CPUs to use during multi processing.
    
    Returns
    -------
    
    """
    # Get all of the schema loci short FASTA files path.
    schema_short_path = os.path.join(schema, 'short')
    schema_loci_short = {loci_path.replace(".fasta", ""): os.path.join(schema_short_path, loci_path) 
                         for loci_path in os.listdir(schema_short_path) 
                         if loci_path.endswith('.fasta')}
    # Create a folder for short translations.
    short_translation_folder = os.path.join(output_path, "short_translation_folder")
    ff.create_directory(short_translation_folder)
    master_loci_short_translation_path = os.path.join(short_translation_folder, "master_short_loci.fasta")
    
    # Translate each short loci and write to master fasta.
    i = 1
    len_short_folder = len(schema_loci_short)
    with open(master_loci_short_translation_path, 'w') as master_fasta:
        
        for loci, loci_short_path in schema_loci_short.items():
            print(f"Translated {i}/{len_short_folder} CDS")
            loci_short_translation_path = os.path.join(short_translation_folder, f"{loci}.fasta")
            i += 1
            fasta_dict = sf.fetch_fasta_dict(loci_short_path, False)
            translation_dict, _, _ = sf.translate_seq_deduplicate(fasta_dict, 
                                                                  loci_short_translation_path,
                                                                  None,
                                                                  constants[5],
                                                                  False)
            
            for loci_id, sequence in translation_dict.items():
                master_fasta.writelines(">"+loci_id+"\n")
                master_fasta.writelines(str(sequence)+"\n")

    # Run BLASTp between new possible loci vs existing loci (all new loci sequences vs loci short FASTA).
    blastp_results_path = os.path.join(output_path, "blastp_results_vs_loci_results")
    ff.create_directory(blastp_results_path)
    bsr_values = {}
    
    # Create query entries.
    for query in outcomes_translations_reps:
        bsr_values[query] = {}
    i = 1
    total = len(outcomes_translations_reps)
    with concurrent.futures.ProcessPoolExecutor(max_workers=cpu) as executor:
        for res in executor.map(bf.run_all_representative_blasts_multiprocessing,
                                outcomes_translations_reps, 
                                repeat('blastp'),
                                repeat(blastp_results_path),
                                repeat(outcomes_translations_reps),
                                repeat(master_loci_short_translation_path)):
            
            filtered_alignments_dict, _, _, _ = af.get_alignments_dict_from_blast_results(res[1], 0, False, False)
            for query, subjects_dict in filtered_alignments_dict.items():
                for subject_id, results in subjects_dict.items():
                    # Score is the largest one between query-subject alignment.
                    # We want the largest score since there may be various matches
                    # alignments, we are interested in knowing overall BSR score
                    # between matches and not for the local alignment.
                    largest_score = 0
                    for entry_id, result in results.items():
                        if result['score'] > largest_score:
                            largest_score = result['score']
                            bsr_values[res[0]].update({subject_id: bf.compute_bsr(result['score'], self_score_dict[query])})
                            
            print(
                f"Running BLASTn for cluster representatives: {res[0]} - {i}/{total}")
            i += 1

def main(schema, output_directory, allelecall_directory, constants, temp_paths, cpu):

    temp_folder = temp_paths[0]
    file_path_cds = temp_paths[1]

    # Verify if the dataset is small, if it is, keep minimum genomes in which
    # specific CDS cluster is present to 5 if not to 1% of the dataset size.
    if not constants[2]:
        count_genomes_path = os.path.join(temp_folder, "1_cds_prediction")
        number_of_genomes = len(os.listdir(count_genomes_path))
        if number_of_genomes <= 20:
            constants[2] = 5
        else:
            constants[2] = round(number_of_genomes * 0.01)
        
    print("Identifying CDS present in the schema...")
    cds_present = os.path.join(temp_folder,"3_cds_preprocess/cds_deduplication/distinct_cds_merged.hashtable")
    # Get dict of CDS and their sequence hashes.
    decoded_sequences_ids = itf.decode_CDS_sequences_ids(cds_present)

    print("Identifying CDS not present in the schema...")
    # Get dict with CDS ids as key and sequence as values.
    not_included_cds = sf.fetch_fasta_dict(file_path_cds, True)
    total_cds = len(not_included_cds)
    print(f"Identified {total_cds} valid CDS not present in the schema.")
    # Filter by size.
    if constants[5]:
        for key, values in list(not_included_cds.items()):
            if len(values) < constants[5]:
                del not_included_cds[key]
        print(f"{len(not_included_cds)}/{total_cds} have size greater or equal to {constants[5]} bp.")
    else:
        constants[5] = 0
        print("No size threshold was applied to the CDS filtering.")

    # Create directories.
    ff.create_directory(output_directory)

    cds_output = os.path.join(output_directory, "1_CDS_processing")
    ff.create_directory(cds_output)
    # This file contains unique CDS.
    cds_not_present_file_path = os.path.join(cds_output, "CDS_not_found.fasta")
    
    # Count the number of CDS present in the schema and write CDS sequence
    # into a FASTA file.
    frequency_cds = {}
    with open(cds_not_present_file_path, 'w+') as cds_not_found:
        for id_, sequence in not_included_cds.items():
            cds_not_found.writelines(">"+id_+"\n")
            cds_not_found.writelines(str(sequence)+"\n")
            
            hashed_seq = sf.seq_to_hash(str(sequence))
            # if CDS sequence is present in the schema count the number of
            # genomes that it is found minus 1 (subtract the first CDS genome).
            if hashed_seq in decoded_sequences_ids:
                frequency_cds[id_] = len(decoded_sequences_ids[hashed_seq]) - 1
            else:
                frequency_cds[id_] = 0
                

    print("Translate and deduplicate unclassified CDS...")
    # Translate the CDS and find unique proteins using hashes, the CDS with
    # the same hash will be added under that hash in protein_hashes.
    cds_not_present_trans_file_path = os.path.join(cds_output, "CDS_not_found_translation.fasta")
    cds_not_present_untrans_file_path = os.path.join(cds_output, "CDS_not_found_untranslated.fasta")
    # Translate and deduplicate protein sequences.
    cds_translation_dict, protein_hashes, _ = sf.translate_seq_deduplicate(not_included_cds,
                                                                           cds_not_present_trans_file_path,
                                                                           cds_not_present_untrans_file_path,
                                                                           constants[5],
                                                                           True)
    # Print additional information about translations and deduplications.
    print(f"{len(cds_translation_dict)}/{len(not_included_cds)} unique protein translations.")

    print("Extracting minimizers for the translated sequences and clustering...")
    # Create variables to store clustering info.
    reps_groups = {}
    clusters = {}
    reps_sequences = {}

    # Sort by size of proteins.
    cds_translation_dict = {k: v for k, v in sorted(cds_translation_dict.items(),
                                                    key=lambda x: len(x[1]),
                                                    reverse=True)}
    # Cluster by minimizers.
    [clusters, reps_sequences, 
     reps_groups, prot_len_dict] = cf.minimizer_clustering(cds_translation_dict,
                                                           5, 5, True, 1, 
                                                           clusters,
                                                           reps_sequences, 
                                                           reps_groups,
                                                           1, constants[3], 
                                                           constants[4],
                                                           True)
    # Print additional information about clustering.
    total_number_clusters = len(clusters)
    print(f"{len(cds_translation_dict)} unique proteins have been clustered into {total_number_clusters} clusters.")
    singleton_cluster = len([cluster for cluster in clusters if len(cluster) == 1])
    print(f"\tOut of those clusters, {singleton_cluster} are singletons")
    print(f"\tOut of those clusters, {total_number_clusters - singleton_cluster} have more than one CDS.")
    
    # Reformat the clusters output, we are interested only in  the ID of cluster members.
    clusters = {cluster_rep: [value[0] for value in values]
                for cluster_rep, values in clusters.items()}
    # For protein hashes get only those that have more than one CDS.
    protein_hashes = {hash_prot: cds_ids for hash_prot, cds_ids in protein_hashes.items()
                      if len(cds_ids) > 1}
    
    # Add also the unique CDS ID that have the same protein as representative.
    for cluster_rep, values in clusters.items():
        for cds_ids in protein_hashes.values():
            # Break since there is only one possible match in protein_hashes.
            if cluster_rep in cds_ids:
                clusters[cluster_rep] + cds_ids[1:]
                break

    print("Filtering clusters...")
    # Get frequency of cluster.
    frequency_cds_cluster = {rep: sum([frequency_cds[entry] for entry in value]) 
                             for rep, value in clusters.items()}
    # Filter cluster by the total sum of CDS that are present in the genomes, based on input value.
    clusters = {rep: cluster_member for rep, cluster_member in clusters.items() 
                if frequency_cds_cluster[rep] >= constants[2]}

    print("Retrieving kmers similiarity and coverage between representatives...")
    reps_kmers_sim = {}
    # Get the representatives protein sequence.
    reps_translation_dict = {rep_id: rep_seq for rep_id, rep_seq in cds_translation_dict.items()
                             if rep_id in clusters}
    # Sort the representative translation dict from largest to smallest.
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
    # Create directories.
    blast_output = os.path.join(output_directory, "2_BLAST_processing")
    ff.create_directory(blast_output)
    
    blastn_output = os.path.join(blast_output, "BLASTn_processing")
    ff.create_directory(blastn_output)
    # Create directory and files path where to write FASTAs.
    representatives_blastn_folder = os.path.join(blastn_output,
                                                "cluster_representatives_fastas")
    ff.create_directory(representatives_blastn_folder)

    representatives_all_fasta_file = os.path.join(representatives_blastn_folder,
                                                  "all_cluster_representatives.fasta")
    # Write files for BLASTn.
    rep_paths_nuc = {}
    # Master file.
    with open(representatives_all_fasta_file, 'w') as all_fasta:
        for cluster_rep_id in clusters:

            all_fasta.writelines(">"+cluster_rep_id+"\n")
            all_fasta.writelines(str(not_included_cds[cluster_rep_id])+"\n")

            rep_fasta_file = os.path.join(representatives_blastn_folder,
                                          f"cluster_rep_{cluster_rep_id}.fasta")
            rep_paths_nuc[cluster_rep_id] = rep_fasta_file
            # Representative file
            with open(rep_fasta_file, 'w') as rep_fasta:
                rep_fasta.writelines(">"+cluster_rep_id+"\n")
                rep_fasta.writelines(str(not_included_cds[cluster_rep_id])+"\n")
    # Create directory
    blastn_results_folder = os.path.join(blastn_output, "blastn_results")
    ff.create_directory(blastn_results_folder)
    # Run BLASTn for all representatives (rep vs all)
    total_reps = len(rep_paths_nuc)
    representative_blast_results = {}
    representative_blast_results_coords_all = {}
    representative_blast_results_coords_pident = {}
    i = 1
    with concurrent.futures.ProcessPoolExecutor(max_workers=cpu) as executor:
        for res in executor.map(bf.run_all_representative_blasts_multiprocessing,
                                clusters,
                                repeat('blastn'),
                                repeat(blastn_results_folder),
                                repeat(rep_paths_nuc),
                                repeat(representatives_all_fasta_file)):

            filtered_alignments_dict, _, alignment_coords_all, alignment_coords_pident = af.get_alignments_dict_from_blast_results(
                res[1], constants[1], True, False)
            # Save the BLASTn results
            representative_blast_results.update(filtered_alignments_dict)
            representative_blast_results_coords_all.update(alignment_coords_all)
            representative_blast_results_coords_pident.update(alignment_coords_pident)

            print(
                f"Running BLASTn for cluster representatives: {res[0]} - {i}/{total_reps}")
            i += 1

    print("Running BLASTp based on BLASTn results matches...")
    # Obtain the list for what BLASTp runs to do, no need to do all vs all as previously.
    # Based on BLASTn results.
    blastp_runs_to_do = {query: itf.flatten_list([[query],[subject[1]['subject']
                                            for subject in subjects.values()]]) 
                         for query, subjects in representative_blast_results.items()}
    
    # Create directories.
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
    # Write the protein FASTA files.
    rep_paths_prot = {}
    rep_matches_prot = {}    
    for query_id, subjects_ids in blastp_runs_to_do.items():
        # First write the representative protein sequence.
        rep_translation_file = os.path.join(representatives_blastp_folder,
                                            f"cluster_rep_translation_{query_id}.fasta")
        rep_paths_prot[query_id] = rep_translation_file
        with open(rep_translation_file, 'w') as trans_fasta:
            trans_fasta.writelines(">"+query_id+"\n")
            trans_fasta.writelines(str(reps_translation_dict[query_id])+"\n")
        # Then write in another file all of the matches for that protein sequence
        # including the representative itself.
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
    self_score_dict = {}
    none_blastp = {}
    # Create query entries
    for query in blastp_runs_to_do:
        bsr_values[query] = {}
        # For self-score
        self_score_dict[query] = {}
    # Run BLASTp between all BLASTn matches (rep vs all its BLASTn matches)  .      
    i = 1
    with concurrent.futures.ProcessPoolExecutor(max_workers=cpu) as executor:
        for res in executor.map(bf.run_all_representative_blasts_multiprocessing,
                                blastp_runs_to_do, 
                                repeat('blastp'),
                                repeat(blastp_results_folder),
                                repeat(rep_paths_prot),
                                rep_matches_prot.values()):
            
            filtered_alignments_dict, self_score, _, _ = af.get_alignments_dict_from_blast_results(res[1], 0, True, True)
            
            # Get IDS of entries that matched with BLASTn but didnt match with BLASTp
            blastn_entries = list(representative_blast_results[res[0]].keys())
            if filtered_alignments_dict:
                blastp_entries = list(filtered_alignments_dict[res[0]].keys())
            else:
                blastp_entries = {}
            if len(blastn_entries) != len(blastp_entries):
                none_blastp[res[0]] = list(set(blastn_entries).symmetric_difference(set(blastp_entries)))

            # Save self-score
            self_score_dict[res[0]] = self_score
            # Since BLAST may find several local aligments choose the largest one to calculate BSR.
            for query, subjects_dict in filtered_alignments_dict.items():
                for subject_id, results in subjects_dict.items():
                    largest_score = 0
                    # Score is the largest one between query-subject alignment.
                    # We want the largest score since there may be various matches
                    # alignments, we are interested in knowing overall BSR score
                    # between matches and not for the local alignment.
                    for entry_id, result in results.items():
                        if result['score'] > largest_score:
                            bsr_values[query].update({subject_id: bf.compute_bsr(result['score'], self_score)})

            print(f"Running BLASTp for cluster representatives matches: {res[0]} - {i}/{total_blasts}")
            i += 1

    add_items_to_results(representative_blast_results, reps_kmers_sim, bsr_values,
                         representative_blast_results_coords_all,
                         representative_blast_results_coords_pident,
                         frequency_cds_cluster)

    print("Filtering BLAST results into classes...")
    results_output = os.path.join(output_directory, "3_Classes_processing")
    ff.create_directory(results_output)
    report_file_path = os.path.join(results_output, "blast_all_matches.tsv")
    
    # Separate results into different classes.
    classes_outcome, drop_list = separate_blastn_results_into_classes(representative_blast_results,
                                                                      constants)
    # Write all of the BLASTn results to a file.
    alignment_dict_to_file(representative_blast_results, report_file_path, 'w')
    
    print("Processing classes...")
    # Process the results_outcome dict and write individual classes to TSV file.
    [cds_to_keep, relationships] = process_classes(representative_blast_results, classes_outcome, drop_list)
    print("Writting classes results to files...")
    write_processed_results_to_file(cds_to_keep, relationships, representative_blast_results,
                                    classes_outcome, results_output)
    
    print("Wrapping up BLAST results...")
    outcomes_translations_reps = wrap_up_blast_results(cds_to_keep, not_included_cds,
                                                       clusters, results_output, constants, cpu)
    print("Reading schema loci short FASTA files...")
    # Create directory
    results_output = os.path.join(output_directory, "4_Schema_processing")
    ff.create_directory(results_output)
    process_schema(schema, outcomes_translations_reps, results_output, 
                   self_score_dict, constants, cpu)