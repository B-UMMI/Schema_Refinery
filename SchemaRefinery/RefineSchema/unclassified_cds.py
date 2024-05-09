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
                       graphical_functions as gf,
                       pandas_functions as pf)
except ModuleNotFoundError:
    from SchemaRefinery.utils import (file_functions as ff,
                                      sequence_functions as sf,
                                      clustering_functions as cf,
                                      blast_functions as bf,
                                      alignments_functions as af,
                                      kmers_functions as kf,
                                      iterable_functions as itf,
                                      graphical_functions as gf,
                                      pandas_functions as pf)

def alignment_dict_to_file(blast_results_dict, file_path, write_type, add_groups_ids = False):
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
    add_groups_ids : bool, optional
        If to add to the header the CDS_group column.

    Returns
    -------
    No return, writes or appends a file at the file_path
    """
    
    header = ['Query\t',
              'Subject\t',
              'Query_length\t',
              'Subject_length\t',
              'Query_start\t',
              'Query_end\t',
              'Subject_start\t',
              'Subject_end\t',
              'Length\t',
              'Score\t',
              'Number_of_gaps\t',
              'Pident\t',
              'Prot_BSR\t',
              'Prot_seq_Kmer_sim\t',
              'Prot_seq_Kmer_cov\t',
              'Cluster_frequency_in_genomes_query_cds\t',
              'Cluster_frequency_in_genomes_subject_cds\t',
              'Global_palign_all_max\t',
              'Global_palign_pident_min\t',
              'Global_palign_pident_max\t',
              'Palign_local_min\t',
              'Class\n']

    if add_groups_ids:
        header[-1] = 'CDS_group\t'
        header.append('Class\n')

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
                         frequency_cds_cluster, loci_ids, add_groups_ids = None):
    """
    Function to add to BLAST results additional information, it adds:
    bsr: value between the two CDS.
    kmers_sim: kmer similarities.
    kmer_cov: kmer coverage.
    frequency_in_genomes_query_cds: how many times that query appears in the schema genomes.
    frequency_in_genomes_subject_cds: how many subject appers in the schema genomes.
    global_palign_all_max: minimum of how much query or subject covers each other.
    global_palign_pident_min: minimum of how much query or subject covers each other, takes into
    account only entries with specific pident value.
    global_palign_pident_max: maximum of how much query or subject covers each other, takes into
    account only entries with specific pident value.
    local_palign_min: minimum of how much query or subject covers each other, takes into
    account only local alignemnt.
                                  
    Parameters
    ----------                                  
    representative_blast_results : 
        Dict that contains representatibes BLAST results.
    reps_kmers_sim : dict
        Dict that contains values for kmer similarities between CDS.
    bsr_values : dict
        Dict that contains BSR values between CDS, may be None.
    representative_blast_results_coords_all : dict
        Dict that contain the coords for all of the entries.
    representative_blast_results_coords_pident : dict
        Dict that contain the coords for all of the entries above certain pident value.
    frequency_cds_cluster : dict
        Dict that contains sum of frequency of that representatives cluster in the
        genomes of the schema.
    loci_ids : bool
        If IDs of loci representatives are included in the frequency_cds_cluster
        they are in this format loci1_x.

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
            if reps_kmers_sim:
                if subject in reps_kmers_sim[query]:
                    sim, cov = reps_kmers_sim[query][subject]
                else:
                    sim = 0
                    cov = 0
            else:
                sim = '-'
                cov = '-'

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
                global_palign_all_max = min(total_length['query'] / result['query_length'],
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
                local_palign_min = min((result['query_end'] - result['query_start'] + 1) / result['query_length'],
                                       (result['subject_end'] - result['subject_start'] + 1) / result['subject_length'])
                # If the alignment is more than 0
                if local_palign_min >= 0:
                    # if IDs of loci representatives are included in the frequency_cds_cluster
                    # they are in this format loci1_x.
                    if loci_ids:
                        query_before = query
                        query = query.split('_')[0]
                    update_dict = {
                        'bsr': bsr,
                        'kmers_sim': sim,
                        'kmers_cov': cov,
                        'frequency_in_genomes_query_cds': frequency_cds_cluster[query],
                        'frequency_in_genomes_subject_cds': frequency_cds_cluster[subject],
                        'global_palign_all_max': global_palign_all_max,
                        'global_palign_pident_min': global_palign_pident_min,
                        'global_palign_pident_max': global_palign_pident_max,
                        'local_palign_min': local_palign_min
                    }
                    # Get the original loci representived ID.
                    if loci_ids:
                        query = query_before
                    representative_blast_results[query][subject][entry_id].update(update_dict)

                    if add_groups_ids:
                        id_ = itf.identify_string_in_dict(subject, add_groups_ids)
                        if not id_:
                            id_ = subject
                        update_dict = {'cds_group': id_}
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
                       '3b',
                       '4']

    # Process results into classes
    for query, rep_blast_result in representative_blast_results.items():
        for id_subject, matches in rep_blast_result.items():
            for id_, blastn_entry in matches.items():
                query_subject_freq = blastn_entry['frequency_in_genomes_query_cds']/blastn_entry['frequency_in_genomes_subject_cds']
                subject_query_freq = blastn_entry['frequency_in_genomes_subject_cds']/blastn_entry['frequency_in_genomes_query_cds']
                
                if blastn_entry['global_palign_all_max'] >= 0.8:
                    # Based on BSR
                    if blastn_entry['bsr'] >= 0.6:
                        add_class_to_dict('1a')
                    # If BSR <0.6 verify if between CDS there ir more than 10x
                    # difference in presence in the schema
                    elif min([query_subject_freq,subject_query_freq]) <= 0.1:
                        add_class_to_dict('1b')
                    # Add two as separate
                    else:
                        add_class_to_dict('1c')
                # Palign < 0.8        
                elif blastn_entry['global_palign_all_max'] > 0.4 and blastn_entry['global_palign_all_max'] <= 0.8:
                    if blastn_entry['pident'] >= constants[1]:
                        # Verify if between CDS there ir more than 10x
                        # difference in presence in the schema
                        if min([query_subject_freq,subject_query_freq]) <= 0.1:
                            add_class_to_dict('2a')
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
                else:
                    add_class_to_dict('4')

    return classes_outcome

def process_classes(representative_blast_results, classes_outcome, all_alleles = None):
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
    cds_joined_cluster : dict, optional
        Dict that contains as keys the IDS of joined clusters or loci, and values
        its elements.

    Returns
    -------
    cds_to_keep : dict     
        Dict of the CDS to keep by each classification.
    relationships : dict
        Dict that contains relationships between various CDS and clusters.
    important_relationships : dict
        Dict that contains as keys the class and values the decisive relatioships
        between loci/CDS.
    drop_list : dict
        Contains the CDS IDs to be removed from further processing for appearing
        fewer time in genomes than their match.
    """
    # Create variables.
    # Variable to add the CDS what will be kept by class.
    cds_to_keep = {}
    # Variable for other ralationships between CDS e.g if two CDS X and Y are 
    # classified as 1a (join) and there are other classfication with them
    # for example Z and X classified as 2b (keep both) then we add to this dict
    # this case for end-user to know.
    other_relationships = {}
    important_relationships = {}
    # List that contains the various CDS with 1a classification to join.
    cluster_to_join = []
    
    # Set of IDs to drop
    drop_list = set()
    processed_case = 0
    main_ids = []
    for class_ in classes_outcome:
        cds_to_keep[class_] = set()
        other_relationships[class_] = []
        important_relationships[class_] = []
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
                # Process all of the cases that have 1a classification.
                # even if they may be in drop_list
                if class_ == '1a':
                    cds_to_keep[class_].update([query, id_subject])
                    cluster_to_join.append([query, id_subject])
                    important_relationships[class_].append([query, id_subject])
                    continue

                retain = []
                # All of the other classifications
                # Find cases that were already processed or to be dropped.
                processed_cases = itf.flatten_list([[c for c in cds] for cds in cds_to_keep.values()])
                processed_cases += drop_list

                if all_alleles:
                    # Change the ID dict if the number of processed cases increases
                    if processed_case != len(processed_cases):
                        processed_case = len(processed_cases)
                        main_ids = set([itf.identify_string_in_dict(group_id, all_alleles) for group_id in processed_cases])
                        main_ids = set(filter(lambda x: x is not None, main_ids))
                        
                    
                    # Don't run the analysis again if one joined CDS or loci already have some results.
                    loci_id = itf.identify_string_in_dict(query, all_alleles)
                    if loci_id in main_ids:
                        continue

                processed = False
                # Get those cases that query and subject were not processed.
                if query not in processed_cases and id_subject not in processed_cases:
                    cds_to_keep[class_].update([query, id_subject])
                    retain = ['r', 'r']
                    processed = True

                # If query was not processed.
                elif query not in processed_cases:
                    cds_to_keep[class_].add(query)
                    retain = ['r']
                    if id_subject in drop_list:
                        retain.append('ad')
                    else:
                        retain.append('ar')
                    processed = True

                # If subject was not processed.
                elif id_subject not in processed_cases:
                    cds_to_keep[class_].add(id_subject)
                    retain = ['r']
                    if query in drop_list:
                        retain.insert(0, 'ad')
                    else:
                        retain.insert(0, 'ar')
                    processed = True

                if class_ in ['1b', '2a'] and processed:
                    blastn_entry = matches[list(matches.keys())[0]]
                    if blastn_entry['frequency_in_genomes_query_cds'] > blastn_entry['frequency_in_genomes_subject_cds']:
                        if id_subject not in cds_to_keep['1a']:
                            if not id_subject in drop_list:
                                retain[1] = 'd'
                                drop_list.add(id_subject)
                            else:
                                retain[1] = 'ad'

                            if id_subject in cds_to_keep[class_]:
                                cds_to_keep[class_].remove(id_subject)
                        else:
                            retain[1] = 'ar'
                    # Remove the query
                    else:
                        if query not in cds_to_keep['1a']:
                            if not query in drop_list:
                                retain[0] = 'd'
                                drop_list.add(query)
                            else:
                                retain[1] = 'ad'
    
                            if query in cds_to_keep[class_]:
                                cds_to_keep[class_].remove(query)
                        else:
                            retain[0] = 'ar'

                if not itf.partially_contains_fragment_of_list([id_subject, query], important_relationships[class_]) and retain:
                    important_relationships[class_].append([query, id_subject, retain])

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
            
    return cds_to_keep, relationships, important_relationships, drop_list

def write_processed_results_to_file(cds_to_keep, relationships, representative_blast_results,
                                    classes_outcome, all_alleles, cds_matched_loci, output_path):
    """
    Write the results from processed_classes into various files.
    
    Parameters
    ----------
    cds_to_keep : dict
        Dict of the CDS to keep by each classification.
    relationships : dict
        Dict that contains relationships between various clusters.
    representative_blast_results : dict
        Dict that contains BLAST results of the representatives with all of the additional
        info.
    classes_outcome : list
        List of list that contains class IDS used in the next function.
    all_alleles : dict
        Dict that contains the loci and joined group as keys and their alleles and 
        elements IDS as values.
        Can be None if no loci are involved.
    cds_matched_loci : dict
        Dict that contains with which loci alleles did CDS and joined group match.
        Can be None if no loci are involved.
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
    add_groups_ids = False
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
            # For the entries that have more than one entry in BLAST results
            # e.g loci1 has alleles loci1_1 and loci1_2
            if all_alleles:
                is_cds = False
                add_groups_ids = True
                cluster_alleles = []
                for entry in cluster:
                    # Skip results in the entries if are retained CDS
                    if entry not in all_alleles or type(entry) == int:
                        if type(entry) == int:
                            cluster = all_alleles[entry]
                        cluster_type = 'CDS_cluster'
                        is_cds = True
                    else:
                        cluster_type = 'loci'
                        cluster_alleles += all_alleles[entry]

                if not is_cds:
                    cluster = cluster_alleles
            else:
                is_cds = True

            if is_cds and cds_matched_loci:
                queries = []
                if type(id_) == int:
                    queries = cds_matched_loci[id_]
                else:
                    for c in cluster:
                        queries += cds_matched_loci[c]
                write_dict = {query : {subject: {id_: entry for id_, entry in entries.items()}
                                       for subject, entries in subjects.items() if subject in cluster}
                              for query, subjects in representative_blast_results.items()
                              if query in queries}
            else:
                write_dict = {query : {subject: {id_: entry for id_, entry in entries.items()}
                                       for subject, entries in subjects.items()}
                              for query, subjects in representative_blast_results.items()
                              if query in cluster}

            report_file_path = os.path.join(blast_by_cluster_output, f"blast_{cluster_type}_{id_}.tsv")
            alignment_dict_to_file(write_dict, report_file_path, 'w', add_groups_ids)
    
    # Create directory.
    partially_contained_relationships = os.path.join(blast_by_cluster_output, "1_partially_contained_relationships")
    ff.create_directory(partially_contained_relationships)
    blast_results_partially_contained = os.path.join(blast_by_cluster_output, "2_blast_results_partially_contained")
    ff.create_directory(blast_results_partially_contained)
    
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
                
            report_file_path = os.path.join(blast_results_partially_contained, f"blast_relationships_to_{cluster_id}.tsv")
            # What write type to use (if file already exists).
            write_type = 'a' if os.path.exists(report_file_path) else 'w'
            # Write BLAST results to file.
            alignment_dict_to_file(write_dict, report_file_path, write_type, add_groups_ids)
            # Path to the relationships report.
            relationships_report_file_path = os.path.join(partially_contained_relationships, f"relationships_to_cluster_{cluster_id}_report.txt")
            
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

    # Create directory.
    joined_cluster_relationships_output = os.path.join(output_path, "blast_results_by_class")
    ff.create_directory(joined_cluster_relationships_output)
    # Write classes to file.
    for class_ in classes_outcome:
        # Fetch all entries with the desired class.
        write_dict = {query : {subject: {id_: entry for id_, entry in entries.items() if entry['class'] == class_}
                               for subject, entries in subjects.items()}
                      for query, subjects in representative_blast_results.items()}

        report_file_path = os.path.join(joined_cluster_relationships_output,
                                        f"blastn_group_{class_}.tsv")
        
        # Write individual class to file.
        alignment_dict_to_file(write_dict,
                               report_file_path,
                               'w',
                               add_groups_ids)
            
def wrap_up_blast_results(cds_to_keep, not_included_cds, clusters, output_path, 
                          constants, drop_list, loci = None, groups_paths_old = None, frequency_cds_cluster = None):
    """
    This function wraps up the results for processing of the unclassified CDSs
    by writing FASTAs files for the possible new loci to include into the schema
    and creates graphs for each results group.
    
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
    constants : list
        Contains the constants to be used in this function.
    drop_list : list
        Contains the CDS IDs to be removed from further processing for appearing
        fewer time in genomes than their match.

    Returns
    -------
    groups_paths_reps : dict
        Dict that contains as Key the ID of each group while the value is the
        path to the FASTA file that contains its nucleotide sequences.
    reps_trans_dict_cds : dict
        Dict that contais the translations of all the CDSs inside the various
        groups.
    """
    # Create directories.
    
    fasta_folder = os.path.join(output_path, "results_fastas")
    ff.create_directory(fasta_folder)
    
    cds_outcome_results_fastas_folder = os.path.join(fasta_folder, "results_group_dna_fastas")
    ff.create_directory(cds_outcome_results_fastas_folder)
    
    cds_outcome_results_reps_fastas_folder = os.path.join(fasta_folder, "results_group_dna_reps_fastas")
    ff.create_directory(cds_outcome_results_reps_fastas_folder)


    # Save CDS that didn't match as seperate variable
    if cds_to_keep.get('Retained_not_matched_by_blastn'):
        Retained_not_matched_by_blastn = cds_to_keep.pop('Retained_not_matched_by_blastn')
    else:
        Retained_not_matched_by_blastn = None

    # Display info about the results obtained from processing the classes.
    # Get the total number of CDS reps considered for classification.
    if loci:
        count_cases = {}
        loci_cases = {}
        cds_cases = {}
        for class_, cds_set in cds_to_keep.items():
            if class_ == '1a':
                loci_cases['1a'] = {}
                cds_cases['1a'] = {}
                for group, cds in cds_set.items():
                    for c in cds:
                        if c in loci:
                            loci_cases['1a'].setdefault(group, []).append(c)
                        else:
                            cds_cases['1a'].setdefault(group, []).append(c)
            else:
                for cds in cds_set:
                    if cds in loci:
                        loci_cases.setdefault(class_, []).append(cds)
                    else:
                        cds_cases.setdefault(class_, []).append(cds)
        for d in drop_list:
            if d in loci:
                loci_cases.setdefault('dropped', []).append(d)
            else:
                cds_cases.setdefault('dropped', []).append(d)
    else:
        count_cases = {}
        for class_, cds_set in cds_to_keep.items():
            if class_ == '1a':
                count_cases[class_] = len(itf.flatten_list(cds_set.values()))
            else:
                count_cases[class_] = len(cds_set)

    if loci:
        for i, printout in enumerate([cds_cases, loci_cases]):
            print(f"Out of {len(groups_paths_old) if i==0 else len(loci)} {'CDSs groups' if i == 0 else 'loci'}:")
            print(f"\t{len(itf.flatten_list(printout.values()))} {'CDSs' if i == 0 else 'loci'}"
                  f" representatives had matches with BLASTn against the {'schema' if i == 0 else 'CDSs'}.")
            for class_, group in printout.items():
                if class_ in ['3b']:
                    print(f"\t\tOut of those groups, {len(group)} {'CDSs groups' if i == 0 else 'loci'}"
                          f" are classified as {class_}"
                          " and not considered further.")
                elif class_ == '1a':
                    print(f"\t\tOut of those groups, {len(itf.flatten_list(printout['1a']))}"
                          f" {'CDSs groups' if i == 0 else 'loci'} are classified as {class_}"
                          f" and are contained in {len(printout['1a'])} joined groups that were retained.")
                elif class_ == 'dropped':
                    print(f"\t\tOut of those {len(group)} {'CDSs groups' if i== 0 else 'loci'}"
                          f" {'were removed from the analysis' if i== 0 else 'are recommended to be replaced with their matched CDS in the schema.'}")
                else:
                    print(f"\t\tOut of those groups, {len(group)} CDS are classified as {class_} and were retained.")
            if i == 0:
                print(f"\t{len(groups_paths_old) - len(itf.flatten_list(printout.values()))}"
                      " didn't have any BLASTn matches so they were retained.\n")
        # Finish writing which cds_to keep
        cds_to_keep
    else:
        # Write info about the classification results.
        print(f"Out of {len(clusters)} clusters:")
        print(f"\t{sum(count_cases.values()) + len(drop_list)} CDS representatives had matches with BLASTn"
              f" which resulted in {len(itf.flatten_list(cds_to_keep.values()))} groups")
        for class_, count in count_cases.items():
            if class_ in ['3b']:
                print(f"\t\tOut of those groups, {count} CDS are classified as {class_}"
                      " and not considered further.")
            elif class_ == '1a':
                print(f"\t\tOut of those groups, {count} CDS are classified as {class_}"
                      f" and are contained in {len(cds_to_keep['1a'])} joined groups"
                      F" that were retained.")
            else:
                print(f"\t\tOut of those groups, {count} CDS are classified as {class_} and were retained.")
        
        print(f"\t\tOut of those {len(drop_list)} CDSs groups were removed from the analysis.")

        if Retained_not_matched_by_blastn:
           print(f"\t{len(Retained_not_matched_by_blastn)} didn't have any BLASTn matches so they were retained.")
          
        cds_to_keep['Retained_not_matched_by_blastn'] = Retained_not_matched_by_blastn
    
    if loci:
        # Write FASTA files for each CDS group to join or retain.
        print("Writting FASTA file for possible new loci...")
        for case_id, cases in enumerate([cds_cases, loci_cases]):
            # Create the directories
            cds_outcome_results = os.path.join(cds_outcome_results_fastas_folder, f"results_{'CDSs' if case_id == 0 else 'loci'}_fastas")
            ff.create_directory(cds_outcome_results)

            id_folder = os.path.join(output_path, "results_IDs")
            ff.create_directory(id_folder)
            id_report_path = os.path.join(id_folder, f"{'CDS_Results' if case_id == 0 else 'Loci_Results'}.tsv")
            # Write dict to TSV
            ff.write_dict_to_tsv(id_report_path, cases)
            
            groups_paths = {}
            groups_paths_reps = {}
            for class_, cds_list in cases.items():
                i = 1
                for cds in cds_list:
                    if class_ == '1a':
                        class_name_cds = f"joined_{i}"
                        i += 1
                    elif class_ == '3a':
                        class_name_cds = f"for_reference_{class_}_{cds}"  
                    elif class_ == '3b':
                        class_name_cds = f"thrash_{class_}_{cds}"
                    elif class_ == 'dropped':
                        class_name_cds = f"dropped_{cds}"
                    else:
                        class_name_cds = f"retained_{class_}_{cds}"
                    
                    file_path = os.path.join(cds_outcome_results, class_name_cds)
                    origin_path = groups_paths_old[cds] if case_id == 0 else loci[cds]
                    ff.copy_file(origin_path, file_path)

        return groups_paths_reps, groups_paths

    else:
        # Write FASTA files for each CDS group to join or retain.
        print("Writting FASTA and additional files for possible new loci...")

        groups_paths = {}
        groups_paths_reps = {}

        for class_, cds_list in cds_to_keep.items():
            for cds in cds_list:
                if class_ == '1a':
                    class_name_cds = f"joined_{cds}"
                elif class_ == '3a':
                    class_name_cds = f"for_reference_{class_}_{cds}"  
                elif class_ == '3b':
                    class_name_cds = f"thrash_{class_}_{cds}"
                else:
                    class_name_cds = f"retained_{class_}_{cds}"
                
                # Create the files paths for all of the FASTAS sequences for the group
                # and only the representatives
                cds_group_fasta_file = os.path.join(cds_outcome_results_fastas_folder, class_name_cds + '.fasta')    
                cds_group_reps_file = os.path.join(cds_outcome_results_reps_fastas_folder, class_name_cds + '.fasta')
                master_file_rep = os.path.join(fasta_folder, 'master_rep_file.fasta')
                # We want to save only retained or joined groups
                if class_ not in ['3a','3b']:
                    groups_paths[cds] = cds_group_fasta_file
                    groups_paths_reps[cds] = cds_group_reps_file
                # to decrease the number of code lines just iterate over string id as a list
                if type(cds) == str:
                    cds = [cds]
                else:
                    cds = cds_to_keep[class_][cds]
                # Write all of the CDS in the group + all of the sequences in the
                # cluster of those representatives.
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
                # Write to master file.
                if not "thrash" in class_name_cds:
                    write_type = 'a' if os.path.exists(master_file_rep) else 'w'
                    with open(master_file_rep, write_type) as fasta_file:
                        for rep_id in cds:
                            fasta_file.writelines(f">{rep_id}\n")
                            fasta_file.writelines(str(not_included_cds[rep_id])+"\n")
    
            # Create directories.
            groups_trans_folder = os.path.join(fasta_folder, "cds_groups_translation")
            ff.create_directory(groups_trans_folder)
            # Translate possible new loci.
            groups_trans = {}
            for key, group_path in groups_paths.items():
                trans_path = os.path.join(groups_trans_folder, os.path.basename(group_path))
                groups_trans[key] = trans_path
                fasta_dict = sf.fetch_fasta_dict(group_path,
                                                 False)
                trans_dict, _, _ = sf.translate_seq_deduplicate(fasta_dict,
                                                                trans_path,
                                                                None,
                                                                constants[5], 
                                                                False,
                                                                False)
        
            # Translate all clusters into their respective cluster file.
            # Create directories
            group_trans_rep_folder = os.path.join(fasta_folder, "cds_groups_translation_reps")
            ff.create_directory(group_trans_rep_folder)
            # Translate possible new loci representatives.
            groups_trans_reps_paths = {}
            reps_trans_dict_cds = {}
            for key, group_path in groups_paths_reps.items():
                trans_path = os.path.join(group_trans_rep_folder, os.path.basename(group_path))
                groups_trans_reps_paths[key] = trans_path
                fasta_dict = sf.fetch_fasta_dict(group_path,
                                                 False)
                trans_dict, _, _ = sf.translate_seq_deduplicate(fasta_dict, 
                                                                trans_path, 
                                                                None,
                                                                constants[5], 
                                                                False,
                                                                False)
                for id_, sequence in trans_dict.items():
                    reps_trans_dict_cds[id_] = sequence
        
        # File to output clusters, their representatives and members.
        cluster_members_output = os.path.join(output_path, 'cluster_members.tsv')
        
        with open(cluster_members_output, 'w') as cluster_members_file:
            # Write header.
            cluster_members_file.write('Cluster_ID\tRepresentatives_IDs'
                                       '\tRep_cluster_members\tFrequency_of_rep\n')
            for class_, cds_list in cds_to_keep.items():
                for cds in cds_list:
                    if class_ == '1a':
                        cluster_members_file.write(str(cds))
                        cds = cds_to_keep[class_][cds]
                    else:
                        cluster_members_file.write(cds)
                        cds = [cds]

                    for rep_id in cds:
                        cluster_members_file.write('\t' + rep_id)
                        cds_ids = [cds_id for cds_id in clusters[rep_id]]
                        for count, cds_id in enumerate(cds_ids):
                            if count == 0:
                                cluster_members_file.write('\t' + cds_id +
                                                           '\t' + str(frequency_cds_cluster[cds_id]) +
                                                           '\n')
                            else:
                                cluster_members_file.write('\t\t' + cds_id + '\n')

        return groups_paths_reps, groups_paths, reps_trans_dict_cds, master_file_rep

def create_graphs(file_path, output_path, filename, other_plots = None):
    """
    Create graphs based on representative_blast_results written inside a TSV file,
    this function creates severall plots related to palign and protein values, with
    the option to create additional plots based on inputs values.
    
    Parameters
    ----------
    file_path : str
        Path to the TSV file.
    output_path : str
        Path to the output directory.
    other_plots : list, optional
        List that contains additional data to create plots.

    Returns
    -------
    Create an HTML file inside the output_path that contains all of the created
    graphs.
    """
    results_output = os.path.join(output_path, "Graph_folder")
    ff.create_directory(results_output)
    
    blast_results_df = ff.import_df_from_file(file_path, '\t')
    
    # Create boxplots
    traces = []
    for column in ['Global_palign_all_max', 'Global_palign_pident_min', 'Global_palign_pident_max', 'Palign_local_min']:
        traces.append(gf.create_violin_plot(y = blast_results_df[column], name = blast_results_df[column].name))
    
    violinplot1 = gf.generate_plot(traces, "Palign Values between BLAST results", "Column", "Palign")
    
    # Create line plot.
    traces = []
    for column in ['Prot_BSR', 'Prot_seq_Kmer_sim', 'Prot_seq_Kmer_cov']:
        traces.append(gf.create_violin_plot(y = blast_results_df[column], name = blast_results_df[column].name))
    
    violinplot2 = gf.generate_plot(traces, "Protein values between BLAST results", "BLAST entries ID", "Columns")
    
    # Create other plots
    extra_plot = []
    if other_plots:
        for plot in other_plots:
            plot_df = pf.dict_to_df(plot[0])
            for column in plot_df.columns.tolist():
                if plot[1] == 'histogram':
                   trace = gf.create_histogram(x = plot_df[column], name = plot_df[column].name)
            
            extra_plot.append(gf.generate_plot(trace, plot[2], plot[3], plot[4]))

    gf.save_plots_to_html([violinplot1, violinplot2] + extra_plot, results_output, filename)

def process_schema(schema, groups_paths, results_output, reps_trans_dict_cds, 
                   cds_to_keep, cds_present, frequency_cds_cluster, allelecall_directory, 
                   master_file_rep, not_included_cds, clusters, constants, cpu):
    """
    This function processes data related to the schema seed, importing, translating
    and BLASTing against the unclassified CDS clusters representatives groups to
    validate them.
    
    Parameters
    ----------
    schema : str
        Path to the schema seed folder.
    groups_trans_reps : dict
        Dict with the paths to the translations of the unclassified CDS clusters.
    results_output : str
        Path were to write the results of this function.
    reps_trans_dict_cds : dict
        Dict that contains the translations for each CDS.
    cds_to_keep : dict     
        Dict of the CDS to keep by each classification.
    cds_present : dict
        Dict that contains the frequency of each CDS in the genomes.
    master_file_rep : str
        Path to the maste file containing retained CDS.
    frequency_cds_cluster : dict
        Contains the frequency of each CDS/loci in the genomes.
    allelecall_directory : str
        Path to the allele call directory.
    not_included_cds : dict
        Dict that contains the unclassified CDS IDs as keys and their
        DNA sequences as values.
    clusters : dict
        Dict that contains the cluster representatives CDS IDs as keys and their
        cluster elements as values.
    constants : list
        Contains the constants to be used in this function.
    cpu : int
        Number of CPUs to use during multi processing.

    Returns
    -------
    representative_blast_results : dict
        Dict that contains BLAST results of the representatives with all of the additional
        info.

    """
    # Get all of the schema loci short FASTA files path.
    schema_short_path = os.path.join(schema, 'short')
    schema_loci_short = {loci_path.replace("_short.fasta", ""): os.path.join(schema_short_path, loci_path) 
                         for loci_path in os.listdir(schema_short_path) 
                         if loci_path.endswith('.fasta')}
    # Create a folder for short translations.
    short_translation_folder = os.path.join(results_output, "short_translation_folder")
    ff.create_directory(short_translation_folder)

    # Find the file in the allele call results that contains the total of each.
    # classification obtained for each loci.
    results_statistics = os.path.join(allelecall_directory, 'loci_summary_stats.tsv')
    # Convert TSV table to dict.
    results_statistics_dict = itf.tsv_to_dict(results_statistics)
    # Add the results for all of the Exact matches to the frequency_cds_cluster dict.
    for key, value in results_statistics_dict.items():
        frequency_cds_cluster.setdefault(key, int(value[0]))
    # Translate each short loci and write to master fasta.
    i = 1
    len_short_folder = len(schema_loci_short)
    all_alleles = {}
    for loci, loci_short_path in schema_loci_short.items():
        print(f"\rTranslated fasta short loci: {i}/{len_short_folder}", end='', flush=True)
        i += 1
        fasta_dict = sf.fetch_fasta_dict(loci_short_path, False)
        
        for loci_id, sequence in fasta_dict.items():
            all_alleles.setdefault(loci, []).append(loci_id)

        loci_short_translation_path = os.path.join(short_translation_folder, f"{loci}.fasta")
        translation_dict, _, _ = sf.translate_seq_deduplicate(fasta_dict, 
                                                              loci_short_translation_path,
                                                              None,
                                                              constants[5],
                                                              False)
        for loci_id, sequence in translation_dict.items():
            reps_trans_dict_cds[loci_id] = sequence

    [representative_blast_results,
     representative_blast_results_coords_all,
     representative_blast_results_coords_pident,
     bsr_values,
     _] = run_blasts(master_file_rep,
                     schema_loci_short,
                     reps_trans_dict_cds,
                     schema_loci_short,
                     results_output,
                     constants,
                     cpu,
                     all_alleles)

    add_items_to_results(representative_blast_results,
                         None,
                         bsr_values,
                         representative_blast_results_coords_all,
                         representative_blast_results_coords_pident,
                         frequency_cds_cluster,
                         True,
                         cds_to_keep['1a'])

    # Add CDS joined clusters to all_alleles IDS
    cds_joined_cluster = cds_to_keep['1a']
    all_alleles.update(cds_joined_cluster)
    # Separate results into different classes.
    classes_outcome = separate_blastn_results_into_classes(representative_blast_results,
                                                           constants)
    
    report_file_path = os.path.join(results_output, "blast_all_matches.tsv")
    # Write all of the BLASTn results to a file.
    alignment_dict_to_file(representative_blast_results, report_file_path, 'w', True)
    
    print("\nProcessing classes...")
    # Process the results_outcome dict and write individual classes to TSV file.
    [cds_to_keep, relationships,
     important_relationships, drop_list] = process_classes(representative_blast_results,
                                                           classes_outcome,
                                                           all_alleles)
    # Replace the alleles entries with their loci ID.
    cds_to_keep = {
        class_: set(
            [entry if not itf.identify_string_in_dict(entry, all_alleles) else itf.identify_string_in_dict(entry, all_alleles) for entry in entries]
        )
        if class_ != '1a' else entries for class_, entries in cds_to_keep.items()
    }

    drop_list = set([entry if not itf.identify_string_in_dict(entry, all_alleles) else itf.identify_string_in_dict(entry, all_alleles) for entry in drop_list])
    # Filter repeated entries
    seen = set()
    for class_, entries in list(cds_to_keep.items()):
        for entry in list(entries):
            if entry not in seen:
                seen.add(entry)
            else:
                cds_to_keep[class_].remove(entry)

    cds_matched_loci = {}
    for class_, entries in list(cds_to_keep.items()):
        for entry in list(entries):
            if entry not in schema_loci_short:
                if type(entry) == int:
                    id_ = entry
                    entry = cds_joined_cluster[entry]
                else:
                    id_ = entry
                    entry = [entry]
                cds_matched_loci.setdefault(id_, set([i[0] for i in itf.flatten_list(important_relationships.values()) if i[1] in entry]))

    print("\nWritting classes results to files...")
    write_processed_results_to_file(cds_to_keep,
                                    relationships,
                                    representative_blast_results,
                                    classes_outcome,
                                    all_alleles,
                                    cds_matched_loci,
                                    results_output)
    
    print("\nWrapping up BLAST results...")
    
    report_main_relationships(important_relationships,
                              representative_blast_results,
                              all_alleles,
                              True,
                              results_output)

    [groups_paths_reps,
     groups_paths] = wrap_up_blast_results(cds_to_keep,
                                           not_included_cds,
                                           all_alleles,
                                           results_output,
                                           constants,
                                           drop_list,
                                           schema_loci_short,
                                           groups_paths,
                                           None)

    return representative_blast_results

def run_blasts(master_fasta_to_blast_against, cds_to_blast, reps_translation_dict,
               rep_paths_nuc, output_dir, constants, cpu, multi_fasta = None):
    """
    This functions runs both BLASTn and Subsequently BLASTp based on results of
    BLASTn.
    
    Parameters
    ----------
    master_fasta_to_blast_against : str
        Path to the master FASTA file that contains all of the nucleotide sequences
        to BLASTn against.
    cds_to_blast : list
        A list that contains all of the ids to BLASTn against master FASTA file.
    reps_translation_dict : dict
        Dict that contains the translations of all the sequences in the master file
        and the CDSs to BLASTn against master file.
    rep_paths_nuc : dict
        Dict that contains the ID of the CDSs to BLASTn against master file
        while the value is the path to the FASTA that contains those CDSs.
    output_dir : str
        Path to write the output of the whole function.
    constants : list
        Contains the constants to be used in this function.
    cpu : int
        Number of CPUs to use during multi processing.
    multi_fasta : dict, optional
        This dict is used as argument in case there are more than one element
        inside each of CDSs FASTA file, since in the initial input of the FASTA
        file may contain more than one CDSs (in case there are more than one
        representatives). BLASTn will perform BLAST with the right file, while
        BLASTp is performed based on results of the BLASTn, since BLASTn result
        are added indiscriminately inside a dict, this multi_fasta dict allows 
        to destinguish which CDS that had matched with BLASTn and to add them inside 
        the common group protein FASTA file to perform BLASTp so the results of 
        BLASTp are more compacted and the results file represent their original 
        input group.
        
    Returns
    -------
    representative_blast_results : 
        Dict that contains representatibes BLAST results.
    representative_blast_results_coords_all : dict
        Dict that contain the coords for all of the entries.
    representative_blast_results_coords_pident : dict
        Dict that contain the coords for all of the entries above certain pident value.
    bsr_values : dict
        Dict that contains BSR values between CDS.
    self_score_dict : dict
        This dict contains the self-score values for all of the CDSs that are
        processed in this function.
    """
    
    print("\nRunning BLASTn between cluster representatives..." if not multi_fasta else
          "\nRunning BLASTn between groups representatives against schema loci short...")
    # BLASTn folder
    blastn_output = os.path.join(output_dir, "BLASTn_processing")
    # Create directory
    blastn_results_folder = os.path.join(blastn_output, "BLASTn_results")
    ff.create_directory(blastn_results_folder)
    # Run BLASTn for all representatives (rep vs all)
    # Calculate max id length for print.
    max_id_length = len(max(cds_to_blast))
    total_reps = len(rep_paths_nuc)
    representative_blast_results = {}
    representative_blast_results_coords_all = {}
    representative_blast_results_coords_pident = {}
    i = 1
    with concurrent.futures.ProcessPoolExecutor(max_workers=cpu) as executor:
        for res in executor.map(bf.run_all_representative_blasts_multiprocessing,
                                cds_to_blast,
                                repeat('blastn'),
                                repeat(blastn_results_folder),
                                repeat(rep_paths_nuc),
                                repeat(master_fasta_to_blast_against)):

            filtered_alignments_dict, _, alignment_coords_all, alignment_coords_pident = af.get_alignments_dict_from_blast_results(
                res[1], constants[1], True, False, True)
            # Save the BLASTn results
            representative_blast_results.update(filtered_alignments_dict)
            representative_blast_results_coords_all.update(alignment_coords_all)
            representative_blast_results_coords_pident.update(alignment_coords_pident)

            print(
                f"\rRunning BLASTn for cluster representatives: {res[0]} - {i}/{total_reps: <{max_id_length}}", 
                end='', flush=True)
            i += 1

    print("\nRunning BLASTp based on BLASTn results matches...")
    # Obtain the list for what BLASTp runs to do, no need to do all vs all as previously.
    # Based on BLASTn results.
    blastp_runs_to_do = {query: itf.flatten_list([[subject[1]['subject']
                                            for subject in subjects.values()]]) 
                         for query, subjects in representative_blast_results.items()}
    
    # Create directories.
    blastp_results = os.path.join(output_dir,
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
    
    blastp_results_ss_folder = os.path.join(blastp_results,
                                            "BLASTp_results_self_score_results")
    ff.create_directory(blastp_results_ss_folder)
    # Write the protein FASTA files.
    rep_paths_prot = {}
    rep_matches_prot = {}
    if multi_fasta:
        blasts_to_run = {}
        seen_entries = {} 
    for query_id, subjects_ids in blastp_runs_to_do.items():
        
        if multi_fasta:
            filename = itf.identify_string_in_dict(query_id, multi_fasta)
            if filename:
                if type(filename) == int:
                    filename = 'joined' + str(filename)
                blasts_to_run.setdefault(filename, set()).update(subjects_ids)
                seen_entries[filename] = set()
            else:
                filename = query_id
                seen_entries[filename] = set()
                blasts_to_run.setdefault(filename, set()).update(subjects_ids)
        else:
            filename = query_id
        # First write the representative protein sequence.
        rep_translation_file = os.path.join(representatives_blastp_folder,
                                            f"cluster_rep_translation_{filename}.fasta")
        
        write_type = 'a' if os.path.exists(rep_translation_file) else 'w'
        
        rep_paths_prot[filename] = rep_translation_file
        with open(rep_translation_file, write_type) as trans_fasta_rep:
            trans_fasta_rep.writelines(">"+query_id+"\n")
            trans_fasta_rep.writelines(str(reps_translation_dict[query_id])+"\n")
        # Then write in another file all of the matches for that protein sequence
        # including the representative itself.
        rep_matches_translation_file = os.path.join(blastn_results_matches_translations,
                                                    f"cluster_matches_translation_{filename}.fasta")
        
        rep_matches_prot[filename] = rep_matches_translation_file
        with open(rep_matches_translation_file, write_type) as trans_fasta:            
            for subject_id in subjects_ids:
                if multi_fasta:
                    if subject_id in seen_entries[filename]:
                        continue

                trans_fasta.writelines(">"+subject_id+"\n")
                trans_fasta.writelines(str(reps_translation_dict[subject_id])+"\n")
                
            if multi_fasta:  
                seen_entries.setdefault(filename, set()).update(subjects_ids)

    # Calculate BSR based on BLASTp.
    bsr_values = {}
    none_blastp = {}
    
    # Create query entries
    for query in blastp_runs_to_do:
        bsr_values[query] = {}
        
    if multi_fasta:
        blastp_runs_to_do = blasts_to_run
    # Total number of runs
    total_blasts = len(blastp_runs_to_do)
    # If there is need to calculate self-score
    print("\nCalculate self-score for the CDSs...")
    self_score_dict = {}
    for query in rep_paths_prot:
        # For self-score
        self_score_dict[query] = {}
    # Calculate self-score
    i = 1
    with concurrent.futures.ProcessPoolExecutor(max_workers=cpu) as executor:
        for res in executor.map(bf.run_self_score_multiprocessing,
                                rep_paths_prot.keys(),
                                repeat('blastp'),
                                rep_paths_prot.values(),
                                repeat(blastp_results_ss_folder)):
            
            _, self_score, _, _ = af.get_alignments_dict_from_blast_results(res[1], 0, True, True, True)
    
            # Save self-score
            self_score_dict[res[0]] = self_score
                            
            print(f"\rRunning BLASTp to calculate self-score for {res[0]: <{max_id_length}}", end='', flush=True)
            i += 1
    # Print newline
    print('\n')  
    
    print("Running BLASTp for representatives..." if not multi_fasta
          else "Running BLASTp for representatives against schema short...")
    # Run BLASTp between all BLASTn matches (rep vs all its BLASTn matches)  .      
    i = 1
    with concurrent.futures.ProcessPoolExecutor(max_workers=cpu) as executor:
        for res in executor.map(bf.run_all_representative_blasts_multiprocessing,
                                blastp_runs_to_do, 
                                repeat('blastp'),
                                repeat(blastp_results_folder),
                                repeat(rep_paths_prot),
                                rep_matches_prot.values()):
            
            filtered_alignments_dict, _, _, _ = af.get_alignments_dict_from_blast_results(res[1], 0, True, False, True)
            
            if not multi_fasta:
                # Get IDS of entries that matched with BLASTn but didnt match with BLASTp
                blastn_entries = list(representative_blast_results[res[0]].keys())
                if filtered_alignments_dict:
                    blastp_entries = list(filtered_alignments_dict[res[0]].keys())
                else:
                    blastp_entries = {}
                if len(blastn_entries) != len(blastp_entries):
                    none_blastp[res[0]] = list(set(blastn_entries).symmetric_difference(set(blastp_entries)))

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
                            largest_score = self_score_dict[res[0]]
                            bsr_values[query].update({subject_id: bf.compute_bsr(result['score'], self_score_dict[res[0]])})
        
            print(f"\rRunning BLASTp for cluster representatives matches: {res[0]} - {i}/{total_blasts: <{max_id_length}}", end='', flush=True)
            i += 1

    return [representative_blast_results, representative_blast_results_coords_all,
            representative_blast_results_coords_pident, bsr_values, self_score_dict]

def report_main_relationships(important_relationships, representative_blast_results,
                              all_alleles, loci, results_output):
    """
    This function reports on the decisive relationships. Based on these relationships the
    decision to which class those CDS or loci were added and if they were ratained
    or removed.

    Parameters
    ----------
    important_relationships : dict
        Dict that contains as keys the class and values the decisive relatioships
        between loci/CDS.
    representative_blast_results : dict
        Dict that contains BLAST results of the representatives with all of the additional
        info.
    all_alleles : dict
        Dict that contains the loci and joined group as keys and their alleles and 
        elements IDS as values.
    loci : bool
        If loci allele IDs are presebt in the important_relationships and representative_blast_results.
    results_output : str
        Path were to write the results of this function.
        
    Returns
    -------
    Creates various TSV files inside the results_output directory for each
    imporant relationship present.
    """
    # Create directories and files
    relationship_output_dir = os.path.join(results_output, "Relationships_results")
    ff.create_directory(relationship_output_dir)
    
    for class_, relationships in important_relationships.items():
        for relationship in relationships:
            query_id = relationship[0]
            # If loci id is present e.g loci1_1.
            if loci:
                id_1 = relationship[0].split('_')[0]
            # CDS ID e.g CDS.
            else:
                # If CDS is part of joined group.
                id_1 = itf.identify_string_in_dict(relationship[0], all_alleles)
                if not id_1:
                    id_1 = relationship[0]
            # If CDS is part of joined group.
            id_2 = itf.identify_string_in_dict(relationship[1], all_alleles)
            subject_id = relationship[1]
            if not id_2:
                id_2 = relationship[1]
            # Get BLAST results by query, subject and class.
            write_dict = {query : {subject: {id_: entry for id_, entry in entries.items() if entry['class'] == class_}
                                   for subject, entries in subjects.items() if subject == subject_id}
                          for query, subjects in representative_blast_results.items() if query == query_id}
            
            # If two CDS were part of the same joined cluster just keep the ID of the cluster.
            file_name = f"{id_1}_vs_{id_2}_{class_}.tsv" if id_1 != id_2 else f"{id_1}_{class_}.tsv"
            
            if class_ != '1a':
                retain = relationship[2]
                
                file_name = retain[0] + '_' + file_name.replace('vs_', 'vs_' + retain[1] + '_')

            report_file_path = os.path.join(relationship_output_dir, file_name)
            # If file already exists just append.
            write_type = 'a' if os.path.exists(report_file_path) else 'w'
            # if add CDS cluster id.
            cds_cluster = True if loci else False
            # Write to file the results.
            alignment_dict_to_file(write_dict, report_file_path, write_type, cds_cluster)

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
    
    # Count CDS size
    cds_size = {}
    for key, sequence in not_included_cds.items():
        cds_size.setdefault(key, len(str(sequence)))

    total_cds = len(not_included_cds)
    print(f"\nIdentified {total_cds} valid CDS not present in the schema.")
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
                

    print("\nTranslate and deduplicate unclassified CDS...")
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
    # Count translation sizes.
    cds_translation_size = {}
    for key, sequence in cds_translation_dict.items():
        cds_translation_size.setdefault(key, len(sequence))

    # Print additional information about translations and deduplications.
    print(f"\n{len(cds_translation_dict)}/{len(not_included_cds)} unique protein translations.")

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
                                                           5,
                                                           5,
                                                           True,
                                                           1, 
                                                           clusters,
                                                           reps_sequences, 
                                                           reps_groups,
                                                           1,
                                                           constants[3], 
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

    print("\nFiltering clusters...")
    # Get frequency of cluster.
    frequency_cds_cluster = {rep: sum([frequency_cds[entry] for entry in value]) 
                             for rep, value in clusters.items()}
    # Filter cluster by the total sum of CDS that are present in the genomes, based on input value.
    clusters = {rep: cluster_member for rep, cluster_member in clusters.items() 
                if frequency_cds_cluster[rep] >= constants[2]}
    print(f"After filtering by CDS frequency in the genomes (>= {constants[2]}),"
          f" out of {total_number_clusters} clusters, {len(clusters)} remained.")

    print("\nRetrieving kmers similiarity and coverage between representatives...")
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
                                                5,
                                                5,
                                                1,
                                                True,
                                                True)
                        )
        
        reps_kmers_sim[cluster_id] = cf.select_representatives(kmers_rep,
                                                               reps_groups,
                                                               0,
                                                               0,
                                                               prot_len_dict,
                                                               cluster_id,
                                                               5)

        reps_kmers_sim[cluster_id] = {match_values[0]: match_values[1:]
                                      for match_values in reps_kmers_sim[cluster_id]}

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
    
    # Run the BLASTn and BLASTp
    [representative_blast_results,
     representative_blast_results_coords_all,
     representative_blast_results_coords_pident,
     bsr_values,
     self_score_dict] = run_blasts(representatives_all_fasta_file,
                                   clusters,
                                   reps_translation_dict,
                                   rep_paths_nuc,
                                   blast_output,
                                   constants,
                                   cpu)
    
    # Add various results to the dict
    add_items_to_results(representative_blast_results,
                         reps_kmers_sim,
                         bsr_values,
                         representative_blast_results_coords_all,
                         representative_blast_results_coords_pident,
                         frequency_cds_cluster,
                         False)

    print("\nFiltering BLAST results into classes...")
    results_output = os.path.join(output_directory, "3_Classes_processing")
    ff.create_directory(results_output)
    report_file_path = os.path.join(results_output, "blast_all_matches.tsv")
    
    # Separate results into different classes.
    classes_outcome = separate_blastn_results_into_classes(representative_blast_results,
                                                           constants)
    # Write all of the BLASTn results to a file.
    alignment_dict_to_file(representative_blast_results, report_file_path, 'w')
    
    print("Processing classes...")
    # Process the results_outcome dict and write individual classes to TSV file.
    [cds_to_keep, relationships, important_relationships, drop_list] = process_classes(representative_blast_results,
                                                                                       classes_outcome)

    report_main_relationships(important_relationships,
                              representative_blast_results,
                              cds_to_keep['1a'],
                              False,
                              results_output)
    
    print("\nAdd remaining cluster that didn't match by BLASTn...")
    # Add cluster not matched by BLASTn
    cds_to_keep['Retained_not_matched_by_blastn'] = set([cluster for cluster in clusters.keys() if cluster not in representative_blast_results.keys()])

    print("\nWritting classes results to files...")
    write_processed_results_to_file(cds_to_keep,
                                    relationships,
                                    representative_blast_results,
                                    classes_outcome,
                                    None,
                                    None,
                                    results_output)
    
    print("\nWrapping up BLAST results...")
    [groups_paths_reps,
     groups_paths,
     reps_trans_dict_cds,
     master_file_rep] = wrap_up_blast_results(cds_to_keep,
                                              not_included_cds,
                                              clusters,
                                              results_output,
                                              constants,
                                              drop_list,
                                              None,
                                              None,
                                              frequency_cds_cluster)
    
    # Add new frequencies in genomes for joined groups
    new_cluster_freq = {}
    for cluster_id, cluster_members in cds_to_keep['1a'].items():
        new_cluster_freq[cluster_id] = 0
        for member in cluster_members:
            new_cluster_freq[cluster_id] += frequency_cds_cluster[member]
        for member in cluster_members:
            frequency_cds_cluster[member] = new_cluster_freq[cluster_id]

    print("Create graphs for the BLAST results...")
    cds_size_dicts = {'IDs': cds_size.keys(),
                      'Size': cds_size.values()}
    cds_translation_size_dicts = {'IDs': cds_size.keys(),
                                  'Size': [int(cds/3) for cds in cds_size.values()]}
    create_graphs(report_file_path,
                  results_output,
                  'All_of_CDS_graphs',
                  [[cds_size_dicts, 'histogram', "Nucleotide Size", 'Size', 'CDS'],
                   [cds_translation_size_dicts, 'histogram','Protein Size' , 'Size', 'CDS']])
    
    for file in ff.get_paths_in_directory(os.path.join(results_output, 'blast_results_by_class')):
        create_graphs(file,
                      results_output,
                      f"graphs_class_{os.path.basename(file).split('_')[-1].replace('.tsv', '')}")

    print("\nReading schema loci short FASTA files...")
    # Create directory
    results_output = os.path.join(output_directory, "4_Schema_processing")
    ff.create_directory(results_output)
    # Create BLASTn_processing directory
    blastn_output = os.path.join(results_output, "BLASTn_processing")
    ff.create_directory(blastn_output)
    # Run Blasts for the found loci against schema short
    representative_blast_results = process_schema(schema,
                                                  groups_paths,
                                                  results_output,
                                                  reps_trans_dict_cds,
                                                  cds_to_keep,
                                                  cds_present,
                                                  frequency_cds_cluster,
                                                  allelecall_directory, 
                                                  master_file_rep,
                                                  not_included_cds,
                                                  clusters,
                                                  constants,
                                                  cpu)