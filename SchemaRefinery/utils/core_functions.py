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
              'Global_palign_all_min\t',
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
    def get_kmer_values(reps_kmers_sim, query, subject):
        """
        Retrieves kmer values for a given query and subject from the reps_kmers_sim dictionary.
        
        Parameters
        ----------
        reps_kmers_sim : dict
            A dictionary containing kmer similarity and coverage values.
        query : str
            The query sequence ID.
        subject : str
            The subject sequence ID.
    
        Returns
        -------
        sim : float or str
            The similarity value if it exists, otherwise returns 0 or '-'.
        cov : float or str
            The coverage value if it exists, otherwise returns 0 or '-'.
        """
        if reps_kmers_sim:
            if subject in reps_kmers_sim[query]:
                sim, cov = reps_kmers_sim[query][subject]
            else:
                sim = 0
                cov = 0
        else:
            sim = '-'
            cov = '-'
        return sim, cov

    def get_bsr_value(bsr_values, query, subject):
        """
        Retrieves BSR value for a given query and subject from the bsr_values dictionary.
        
        Parameters
        ----------
        bsr_values : dict
            A dictionary containing BSR values.
        query : str
            The query sequence ID.
        subject : str
            The subject sequence ID.
            
        Returns
        -------
        bsr : float
            The BSR value if it exists, otherwise returns 0. If BSR value is greater than 1, it is rounded to the nearest integer.
        """
        bsr = bsr_values[query].get(subject, 0)
        if bsr > 1.0:
            bsr = float(round(bsr))
        return bsr

    def calculate_total_length(representative_blast_results_coords, query, subject):
        """
        Calculates total length for a given query and subject.
        
        Parameters
        ----------
        representative_blast_results_coords : dict
            A dictionary containing BLAST results coordinates.
        query : str
            The query sequence ID.
        subject : str
            The subject sequence ID.
            
        Returns
        -------
        total_length : dict
            A dictionary with the total length of each reference sequence.
        """
        total_length = {}
        for ref, intervals in representative_blast_results_coords[query][subject].items():
            if intervals:
                sorted_intervals = sorted(intervals, key=lambda x: x[0])
                length = sum(interval[1] - interval[0] + 1 for interval in af.merge_intervals(sorted_intervals))
                total_length[ref] = length
            else:
                total_length[ref] = 0
        return total_length

    def calculate_global_palign(total_length, result):
        """
        Calculates global palign for a given total length and result.
        
        Parameters
        ----------
        total_length : dict
            A dictionary containing the total length of each reference sequence.
        result : dict
            A dictionary containing the result of a BLAST search.
            
        Returns
        -------
        global_palign_min : float
            The minimum global palign value.
        global_palign_max : float
            The maximum global palign value.
        """
        global_palign_min = min(total_length['query'] / result['query_length'],
                                total_length['subject'] / result['subject_length'])
        global_palign_max = max(total_length['query'] / result['query_length'],
                                total_length['subject'] / result['subject_length'])
        return global_palign_min, global_palign_max

    def calculate_local_palign(result):
        """
        Calculates local palign for a given result.
        
        Parameters
        ----------
        result : dict
            A dictionary containing the result of a BLAST search.
            
        Returns
        -------
        local_palign_min : float
            The minimum local palign value.
        """
        local_palign_min = min((result['query_end'] - result['query_start'] + 1) / result['query_length'],
                            (result['subject_end'] - result['subject_start'] + 1) / result['subject_length'])
        return local_palign_min

    def update_results(representative_blast_results, query, subject, entry_id, bsr, sim, cov, frequency_cds_cluster, global_palign_all_min, global_palign_all_max, global_palign_pident_min, global_palign_pident_max, local_palign_min, loci_ids, add_groups_ids):
        """
        Updates results for a given query and subject.
        
        Parameters
        ----------
        representative_blast_results : dict
            A dictionary containing BLAST results.
        query : str
            The query sequence ID.
        subject : str
            The subject sequence ID.
        entry_id : str
            The ID of the entry to update.
        bsr, sim, cov : float
            The BSR, similarity, and coverage values.
        frequency_cds_cluster : dict
            A dictionary containing the frequency of CDS clusters.
        global_palign_all_min, global_palign_all_max, global_palign_pident_min, global_palign_pident_max, local_palign_min : float
            The minimum and maximum global palign values, and the minimum local palign value.
        loci_ids : list
            A list of loci IDs.
        add_groups_ids : list
            A list of additional group IDs.
        
        Returns
        -------
        No returns, modifies the representative_blast_results dict inside the parent function.
        """
        if loci_ids:
            query_before = query
            query = query.split('_')[0]
        update_dict = {
            'bsr': bsr,
            'kmers_sim': sim,
            'kmers_cov': cov,
            'frequency_in_genomes_query_cds': frequency_cds_cluster[query],
            'frequency_in_genomes_subject_cds': frequency_cds_cluster[subject],
            'global_palign_all_min' : global_palign_all_min,
            'global_palign_all_max': global_palign_all_max,
            'global_palign_pident_min': global_palign_pident_min,
            'global_palign_pident_max': global_palign_pident_max,
            'local_palign_min': local_palign_min
        }
        if loci_ids:
            query = query_before
        representative_blast_results[query][subject][entry_id].update(update_dict)

        if add_groups_ids:
            id_ = itf.identify_string_in_dict(subject, add_groups_ids)
            if not id_:
                id_ = subject
            update_dict = {'cds_group': id_}
            representative_blast_results[query][subject][entry_id].update(update_dict)

    def remove_results(representative_blast_results, query, subject, entry_id):
        """
        Removes results for a given query and subject.
        
        Parameters
        ----------
        representative_blast_results : dict
            A dictionary containing BLAST results.
        query : str
            The query sequence ID.
        subject : str
            The subject sequence ID.
        entry_id : str
            The ID of the entry to remove.
        
        Returns
        -------
        No returns, modifies the representative_blast_results dict inside the parent function.
        """
        del representative_blast_results[query][subject][entry_id]

    def clean_up_results(representative_blast_results, query, subject):
        """
        Cleans up results for a given query and subject by removing empty entries.
        
        Parameters
        ----------
        representative_blast_results : dict
            A dictionary containing BLAST results.
        query : str
            The query sequence ID.
        subject : str
            The subject sequence ID.
        
        Returns
        -------
        No returns, modifies the representative_blast_results dict inside the parent function.
        """
        if not representative_blast_results[query][subject]:
            del representative_blast_results[query][subject]
        if not representative_blast_results[query]:
            del representative_blast_results[query]

    # Iterate over the representative_blast_results dictionary
    for query, subjects_dict in list(representative_blast_results.items()):
        for subject, blastn_results in list(subjects_dict.items()):
            sim, cov = get_kmer_values(reps_kmers_sim, query, subject)
            bsr = get_bsr_value(bsr_values, query, subject)

            # Iterate over the blastn_results dictionary
            for entry_id, result in list(blastn_results.items()):
                total_length = calculate_total_length(representative_blast_results_coords_all, query, subject)
                global_palign_all_min, global_palign_all_max = calculate_global_palign(total_length, result)

                total_length = calculate_total_length(representative_blast_results_coords_pident, query, subject)
                global_palign_pident_min, global_palign_pident_max = calculate_global_palign(total_length, result)

                local_palign_min = calculate_local_palign(result)

                if local_palign_min >= 0:
                    update_results(representative_blast_results, query, subject, entry_id, bsr, sim, cov, frequency_cds_cluster, global_palign_all_min, global_palign_all_max, global_palign_pident_min, global_palign_pident_max, local_palign_min, loci_ids, add_groups_ids)
                else:
                    remove_results(representative_blast_results, query, subject, entry_id)

            clean_up_results(representative_blast_results, query, subject)

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

    # Define classes based on priority
    classes_outcome = ['1a', '1b', '2a', '3a', '2b', '1c', '3b', '4a', '4b', '4c','5']

    # Loop through the representative BLAST results
    for query, rep_blast_result in representative_blast_results.items():
        for id_subject, matches in rep_blast_result.items():
            for id_, blastn_entry in matches.items():
                # Calculate the frequency ratio
                freq_ratio = min(blastn_entry['frequency_in_genomes_query_cds']/blastn_entry['frequency_in_genomes_subject_cds'],
                                blastn_entry['frequency_in_genomes_subject_cds']/blastn_entry['frequency_in_genomes_query_cds'])
                
                # Classify based on global_palign_all_min and bsr
                if blastn_entry['global_palign_all_min'] >= 0.8:
                    if blastn_entry['bsr'] >= 0.6:
                        # Add to class '1a' if bsr is greater than or equal to 0.6
                        add_class_to_dict('1a')
                    elif freq_ratio <= 0.1:
                        # Add to class '1b' if frequency ratio is less than or equal to 0.1
                        add_class_to_dict('1b')
                    else:
                        # Add to class '1c' if none of the above conditions are met
                        add_class_to_dict('1c')
                elif 0.4 <= blastn_entry['global_palign_all_min'] < 0.8:
                    if blastn_entry['pident'] >= constants[1]:
                        if blastn_entry['global_palign_pident_max'] >= 0.8:
                            # Add to class '2a' or '2b' based on frequency ratio
                            add_class_to_dict('2a' if freq_ratio <= 0.1 else '2b')
                        else:
                            # Add to class '3a' or '3b' based on frequency ratio
                            add_class_to_dict('3a' if freq_ratio <= 0.1 else '3b')
                    else:
                        if blastn_entry['global_palign_pident_max'] >= 0.8:
                            # Add to class '4a' or '4b' based on frequency ratio
                            add_class_to_dict('4a' if freq_ratio <= 0.1 else '4b')
                        else:
                            # Add to class '4c' if none of the above conditions are met
                            add_class_to_dict('4c')
                else:
                    # Add to class '5' for everything that is unrelated
                    add_class_to_dict('5')

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
    important_relationships : dict
        Dict that contains as keys the class and values the decisive relatioships
        between loci/CDS.
    drop_list : dict
        Contains the CDS IDs to be removed from further processing for appearing
        fewer time in genomes than their match.
    """
    # Initialize variables
    cds_to_keep = {}
    important_relationships = {}
    cluster_to_join = []
    drop_list = set()
    processed_case = 0
    main_ids = []

    # Loop over each class
    for class_ in classes_outcome:
        cds_to_keep[class_] = set()
        important_relationships[class_] = []

        # Get all entries with the desired class and remove empty entries
        class_dict = {query: {subject: {id_: entry for id_, entry in entries.items() if entry['class'] == class_}
                              for subject, entries in subjects.items()}
                      for query, subjects in representative_blast_results.items()}
        itf.remove_empty_dicts_recursive(class_dict)

        # Process the CDS to find what CDS to retain while also adding the relationships between different CDS
        for query, rep_blast_result in class_dict.items():
            for id_subject, matches in rep_blast_result.items():
                # Process all of the cases that have 1a classification
                if class_ == '1a':
                    cds_to_keep[class_].update([query, id_subject])
                    cluster_to_join.append([query, id_subject])
                    important_relationships[class_].append([query, id_subject])
                    continue

                # Initialize retain list
                retain = []

                # Find cases that were already processed or to be dropped
                processed_cases = itf.flatten_list([[c for c in cds] for cds in cds_to_keep.values()])
                processed_cases += drop_list

                # Change the ID dict if the number of processed cases increases
                if all_alleles and processed_case != len(processed_cases):
                    processed_case = len(processed_cases)
                    main_ids = set([itf.identify_string_in_dict(group_id, all_alleles) for group_id in processed_cases])
                    main_ids = set(filter(lambda x: x is not None, main_ids))

                # Don't run the analysis again if one joined CDS or loci already have some results
                if all_alleles and itf.identify_string_in_dict(query, all_alleles) in main_ids:
                    continue

                # Process cases where neither query nor subject were processed
                if query not in processed_cases and id_subject not in processed_cases:
                    cds_to_keep[class_].update([query, id_subject])
                    retain = ['r', 'r']

                # Process cases where only query was not processed
                elif query not in processed_cases:
                    cds_to_keep[class_].add(query)
                    retain = ['r', 'ad' if id_subject in drop_list else 'ar']

                # Process cases where only subject was not processed
                elif id_subject not in processed_cases:
                    cds_to_keep[class_].add(id_subject)
                    retain = ['ad' if query in drop_list else 'ar', 'r']

                if class_ in ['1b', '2a', '3a'] and retain:
                    blastn_entry = matches[list(matches.keys())[0]]
                    is_frequency_greater = blastn_entry['frequency_in_genomes_query_cds'] >= blastn_entry['frequency_in_genomes_subject_cds']
                    id_or_query = id_subject if is_frequency_greater else query
                    retain_index = 1 if is_frequency_greater else 0
                    retain_value = 'ar' if id_or_query in cds_to_keep['1a'] else 'ad' if id_or_query in drop_list else 'd'
                    
                    retain[retain_index] = retain_value
                    if retain_value != 'ar':
                        drop_list.add(id_or_query)
                        cds_to_keep[class_].discard(id_or_query)

                if not itf.partially_contains_fragment_of_list([id_subject, query], important_relationships[class_]) and retain:
                    important_relationships[class_].append([query, id_subject, retain])

    # Create the joined cluster by joining by IDs.
    cds_to_keep['1a'] = {i+1: join for i, join in enumerate(cf.cluster_by_ids(cluster_to_join))}

    return cds_to_keep, important_relationships, drop_list

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
    def print_classification_results(class_, count, printout, i):
        """
        Prints the classification results based on the class type.

        Parameters
        ----------
        class_ : str
            The class type.
        count : int
            The count of groups.
        printout : dict
            The dictionary containing printout information.
        i : int
            An index used to determine the printout message.

        Returns
        -------
        None, prints in stdout
        """
        if count > 0:
            if class_ in ['2b', '4b']:
                print(f"\t\tOut of those groups, {count} {'CDSs' if i == 0 else 'loci'} are classified as {class_} and were retained"
                    " but it is recommended to verify them as they may be contained or contain partially inside"
                    " their BLAST match.")
            elif class_ == '1a':
                print(f"\t\tOut of those groups, {count} {'CDSs groups' if i == 0 else 'loci'} are classified as {class_}"
                    f" and are contained in {len(printout['1a'])} joined groups that were retained.")
            elif class_ == 'dropped':
                print(f"\t\tOut of those {count} {'CDSs groups' if i== 0 else 'loci'}"
                    f" {'were removed from the analysis' if i== 0 else 'are recommended to be replaced with their matched CDS in the schema.'}")
            else:
                print(f"\t\tOut of those groups, {count} {'CDSs' if i == 0 else 'loci'} are classified as {class_} and were retained.")

    def create_directory_and_write_dict(cds_outcome_results_fastas_folder, output_path, case_id, cases):
        """
        Create directories and write dict to TSV.

        Parameters
        ----------
        cds_outcome_results_fastas_folder : str
            The path to the folder where the results will be stored.
        output_path : str
            The path where the output will be written.
        case_id : int
            The ID of the case.
        cases : dict
            The dictionary containing the cases.

        Returns
        -------
        cds_outcome_results : str
            The path to the results folder.
        """
        cds_outcome_results = os.path.join(cds_outcome_results_fastas_folder, f"results_{'CDSs' if case_id == 0 else 'loci'}_fastas")
        ff.create_directory(cds_outcome_results)

        id_folder = os.path.join(output_path, "results_IDs")
        ff.create_directory(id_folder)
        id_report_path = os.path.join(id_folder, f"{'CDS_Results' if case_id == 0 else 'Loci_Results'}.tsv")
        ff.write_dict_to_tsv(id_report_path, cases)

        return cds_outcome_results

    def copy_fasta(class_, cds_list, case_id, cds_outcome_results, groups_paths_old, loci):
        """
        Process each class and CDS list in cases.

        Parameters
        ----------
        class_ : str
            The class type.
        cds_list : list
            The list of CDSs.
        case_id : int
            The ID of the case.
        cds_outcome_results : str
            The path to the results folder.
        groups_paths_old : dict
            The dictionary containing the old paths.
        loci : dict
            The dictionary containing the loci.

        Returns
        -------
        None, copies file from origin to destination.
        """
        i = 1
        for cds in cds_list:
            if class_ == '1a':
                class_name_cds = f"joined_{i}"
                i += 1
            elif class_ == 'dropped':
                class_name_cds = f"dropped_{cds}"
            else:
                class_name_cds = f"retained_{class_}_{cds}"
            
            file_path = os.path.join(cds_outcome_results, class_name_cds)
            origin_path = groups_paths_old.pop(cds) if case_id == 0 else loci[cds]
            ff.copy_file(origin_path, file_path)

    def write_fasta_to_keep(class_, cds_list, cds_outcome_results_fastas_folder, cds_outcome_results_reps_fastas_folder, fasta_folder, groups_paths, groups_paths_reps, not_included_cds, clusters):
        """
        Process each class and CDS list in cds_to_keep.

        Parameters
        ----------
        class_ : str
            The class type.
        cds_list : list
            The list of CDSs.
        cds_outcome_results_fastas_folder : str
            The path to the results folder.
        cds_outcome_results_reps_fastas_folder : str
            The path to the folder where the representative results will be stored.
        fasta_folder : str
            The path to the folder where the fasta files are stored.
        groups_paths : dict
            The dictionary containing the paths to the groups.
        groups_paths_reps : dict
            The dictionary containing the paths to the representative groups.
        not_included_cds : dict
            The dictionary containing the CDSs that were not included.
        clusters : dict
            The dictionary containing the clusters.

        Returns
        -------
        None, writtes FASTA files.
        """
        for cds in cds_list:
            if class_ == '1a':
                class_name_cds = f"joined_{cds}"
            elif class_ == 'Retained_not_matched_by_blastn':
                class_name_cds = f"retained_not_matched_by_blastn_{cds}"
            else:
                class_name_cds = f"retained_{class_}_{cds}"
            
            cds_group_fasta_file = os.path.join(cds_outcome_results_fastas_folder, class_name_cds + '.fasta')    
            cds_group_reps_file = os.path.join(cds_outcome_results_reps_fastas_folder, class_name_cds + '.fasta')
            master_file_rep = os.path.join(fasta_folder, 'master_rep_file.fasta')
            groups_paths[cds] = cds_group_fasta_file
            groups_paths_reps[cds] = cds_group_reps_file
            if type(cds) == str:
                cds = [cds]
            else:
                cds = cds_to_keep[class_][cds]
            with open(cds_group_fasta_file, 'w') as fasta_file:
                for rep_id in cds:
                    cds_ids = [cds_id for cds_id in clusters[rep_id]]
                    for cds_id in cds_ids:
                        fasta_file.writelines(f">{cds_id}\n")
                        fasta_file.writelines(str(not_included_cds[cds_id])+"\n")
            with open(cds_group_reps_file, 'w') as fasta_file:
                for rep_id in cds:
                    fasta_file.writelines(f">{rep_id}\n")
                    fasta_file.writelines(str(not_included_cds[rep_id])+"\n")
            write_type = 'a' if os.path.exists(master_file_rep) else 'w'
            with open(master_file_rep, write_type) as fasta_file:
                for rep_id in cds:
                    fasta_file.writelines(f">{rep_id}\n")
                    fasta_file.writelines(str(not_included_cds[rep_id])+"\n")

    def translate_possible_new_loci(fasta_folder, groups_paths, groups_paths_reps, constants):
        """
        Translate possible new loci and writes to master file.

        Parameters
        ----------
        fasta_folder : str
            The path to the folder where the fasta files are stored.
        groups_paths : dict
            The dictionary containing the paths to the groups.
        groups_paths_reps : dict
            The dictionary containing the paths to the representative groups.
        constants : list
            The list of constants.

        Returns
        -------
        reps_trans_dict_cds :dict
            The dictionary containing the translated sequences.
        """
        groups_trans_folder = os.path.join(fasta_folder, "cds_groups_translation")
        ff.create_directory(groups_trans_folder)
        groups_trans = {}
        for key, group_path in groups_paths.items():
            trans_path = os.path.join(groups_trans_folder, os.path.basename(group_path))
            groups_trans[key] = trans_path
            fasta_dict = sf.fetch_fasta_dict(group_path, False)
            trans_dict, _, _ = sf.translate_seq_deduplicate(fasta_dict, trans_path, None, constants[5], False, False)
        
        group_trans_rep_folder = os.path.join(fasta_folder, "cds_groups_translation_reps")
        ff.create_directory(group_trans_rep_folder)
        groups_trans_reps_paths = {}
        reps_trans_dict_cds = {}
        for key, group_path in groups_paths_reps.items():
            trans_path = os.path.join(group_trans_rep_folder, os.path.basename(group_path))
            groups_trans_reps_paths[key] = trans_path
            fasta_dict = sf.fetch_fasta_dict(group_path, False)
            trans_dict, _, _ = sf.translate_seq_deduplicate(fasta_dict, trans_path, None, constants[5], False, False)
            for id_, sequence in trans_dict.items():
                reps_trans_dict_cds[id_] = sequence

        return reps_trans_dict_cds

    def write_cluster_members_to_file(output_path, cds_to_keep, clusters, frequency_cds_cluster):
        """
        Write cluster members to file.

        Parameters
        ----------
        output_path : str
            The path where the output will be written.
        cds_to_keep : dict
            The dictionary containing the CDSs to keep.
        clusters : dict
            The dictionary containing the clusters.
        frequency_cds_cluster : dict
            The dictionary containing the frequency of each CDS in the cluster.

        Returns
        -------
        None, writes to file.
        """
        cluster_members_output = os.path.join(output_path, 'cluster_members.tsv')
        with open(cluster_members_output, 'w') as cluster_members_file:
            cluster_members_file.write('Cluster_ID\tRepresentatives_IDs\tRep_cluster_members\tFrequency_of_rep\n')
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
                                cluster_members_file.write('\t' + cds_id + '\t' + str(frequency_cds_cluster[cds_id]) + '\n')
                            else:
                                cluster_members_file.write('\t\t' + cds_id + '\n')
    # Create directories.
    
    fasta_folder = os.path.join(output_path, "results_fastas")
    ff.create_directory(fasta_folder)
    
    cds_outcome_results_fastas_folder = os.path.join(fasta_folder, "results_group_dna_fastas")
    ff.create_directory(cds_outcome_results_fastas_folder)
    if not loci:
        cds_outcome_results_reps_fastas_folder = os.path.join(fasta_folder, "results_group_dna_reps_fastas")
        ff.create_directory(cds_outcome_results_reps_fastas_folder)


    # If 'Retained_not_matched_by_blastn' exists in cds_to_keep, remove it and store it separately
    Retained_not_matched_by_blastn = cds_to_keep.pop('Retained_not_matched_by_blastn', None)

    # Display info about the results obtained from processing the classes.
    # Get the total number of CDS reps considered for classification.
    count_cases = {}
    loci_cases = {}
    cds_cases = {}
    if loci:
        # Iterate over classes and their associated CDS sets
        for class_, cds_set in cds_to_keep.items():
            # Initialize dictionaries for class '1a'
            if class_ == '1a':
                loci_cases['1a'] = {}
                cds_cases['1a'] = {}
                # Iterate over groups and their associated CDS in class '1a'
                for group, cds in cds_set.items():
                    # Separate CDS into those in loci and those not in loci
                    loci_cases['1a'][group] = [c for c in cds if c in loci]
                    cds_cases['1a'][group] = [c for c in cds if c not in loci]
            else:
                # For other classes, separate CDS into those in loci and those not in loci
                loci_cases[class_] = [cds for cds in cds_set if cds in loci]
                cds_cases[class_] = [cds for cds in cds_set if cds not in loci]

        # Process drop_list in the same way as above
        loci_cases['dropped'] = [d for d in drop_list if d in loci]
        cds_cases['dropped'] = [d for d in drop_list if d not in loci]

    else:
        for class_, cds_set in cds_to_keep.items():
            if class_ == '1a':
                count_cases[class_] = len(itf.flatten_list(cds_set.values()))
            else:
                count_cases[class_] = len(cds_set)

    # Check if loci is not empty
    if loci:
        for i, printout in enumerate([cds_cases, loci_cases]):
            print(f"Out of {len(groups_paths_old) if i==0 else len(loci)} {'CDSs groups' if i == 0 else 'loci'}:")
            print(f"\t{len(itf.flatten_list(printout.values()))} {'CDSs' if i == 0 else 'loci'}"
                f" representatives had matches with BLASTn against the {'schema' if i == 0 else 'CDSs'}.")

            # Print the classification results
            for class_, group in printout.items():
                print_classification_results(class_, len(group), printout, i)

            if i == 0:
                print(f"\t{len(groups_paths_old) - len(itf.flatten_list(printout.values()))}"
                    " didn't have any BLASTn matches so they were retained.\n")
    else:
        # Write info about the classification results.
        print(f"Out of {len(clusters)} clusters:")
        print(f"\t{sum(count_cases.values()) + len(drop_list)} CDS representatives had matches with BLASTn"
            f" which resulted in {len(itf.flatten_list(cds_to_keep.values()))} groups")

        # Print the classification results
        for class_, count in count_cases.items():
            print_classification_results(class_, count, cds_to_keep, 0)

        print(f"\t\tOut of those {len(drop_list)} CDSs groups were removed from the analysis.")

        if Retained_not_matched_by_blastn:
            print(f"\t{len(Retained_not_matched_by_blastn)} didn't have any BLASTn matches so they were retained.")
            
            cds_to_keep['Retained_not_matched_by_blastn'] = Retained_not_matched_by_blastn

    # Initialize dictionaries to store paths
    groups_paths = {}
    groups_paths_reps = {}
    # Check if loci is not None.
    if loci:
        print("Writing FASTA file for possible new loci...")
        for case_id, cases in enumerate([cds_cases, loci_cases]):
            # Create directories and write dict to TSV
            cds_outcome_results = create_directory_and_write_dict(cds_outcome_results_fastas_folder, output_path, case_id, cases)

            # Process each class and CDS list in cases
            for class_, cds_list in cases.items():
                copy_fasta(class_, cds_list, case_id, cds_outcome_results, groups_paths_old, loci)
            # Copy CDS that didnt match
            for cds, path in groups_paths_old.items():
                cds_name = f"retained_not_matched_by_blastn_{cds}"
                file_path = os.path.join(cds_outcome_results, cds_name)
                ff.copy_file(path, file_path)

        master_file_rep = None
        reps_trans_dict_cds = None
    else:
        print("Writing FASTA and additional files for possible new loci...")

        # Process each class and CDS list in cds_to_keep
        for class_, cds_list in cds_to_keep.items():
            write_fasta_to_keep(class_, cds_list, cds_outcome_results_fastas_folder, cds_outcome_results_reps_fastas_folder, fasta_folder, groups_paths, groups_paths_reps, not_included_cds, clusters)

        # Translate possible new loci and write to master file
        reps_trans_dict_cds = translate_possible_new_loci(fasta_folder, groups_paths, groups_paths_reps, constants)

        # Write cluster members to file
        write_cluster_members_to_file(output_path, cds_to_keep, clusters, frequency_cds_cluster)

        master_file_rep = os.path.join(fasta_folder, 'master_rep_file.fasta')
    return groups_paths, reps_trans_dict_cds, master_file_rep

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
