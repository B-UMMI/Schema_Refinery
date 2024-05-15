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
    classes_outcome = ['1a', '1b', '2a', '2b', '3a', '1c', '3b', '4a', '4b', '4c', '5a', '5b', '5c']

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
                    if blastn_entry['global_palign_pident_max'] >= 0.8:
                        # Add to class '5a' or '5b' based on frequency ratio
                        add_class_to_dict('5a' if freq_ratio <= 0.1 else '5b')
                    else:
                        # Add to class '5c' if none of the above conditions are met
                        add_class_to_dict('5c')

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

                if class_ in ['1b', '2a', '3a', '5a'] and retain:
                    blastn_entry = matches[list(matches.keys())[0]]
                    is_frequency_greater = blastn_entry['frequency_in_genomes_query_cds'] > blastn_entry['frequency_in_genomes_subject_cds']
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

def write_processed_results_to_file(cds_to_keep, representative_blast_results,
                                    classes_outcome, all_alleles, cds_matched_loci, output_path):
    """
    Write the results from processed_classes into various files.
    
    Parameters
    ----------
    cds_to_keep : dict
        Dict of the CDS to keep by each classification.
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
    def process_clusters(cds_to_keep, representative_blast_results, all_alleles, cds_matched_loci, output_path):
        """
        Process and write cluster results.

        Parameters
        ----------
        cds_to_keep : dict
            Dictionary of classes and their corresponding CDS.
        representative_blast_results : dict
            Dictionary of representative blast results.
        all_alleles : dict
            Dict that contains the IDs as key and of all alleles related to that ID as values.
        cds_matched_loci : dict
            Dictionary of CDS matched loci.
        output_path : str
            Path to the output directory.

        Returns
        -------
        add_groups_ids : bool
            True if additional group IDs are present, False otherwise.
        """
        # Loop over each class and its corresponding CDS
        for class_, cds in cds_to_keep.items():
            # Loop over each cluster in the CDS
            for id_, cluster in enumerate(cds, 1):
                # Process the cluster and get the necessary details
                id_, cluster, cluster_type, is_cds, add_groups_ids = process_cluster(class_, id_, cluster, all_alleles, cds)
                # Generate a dictionary to be written to the file
                write_dict = generate_write_dict(id_, cluster, is_cds, cds_matched_loci, representative_blast_results)
                # Define the path of the report file
                report_file_path = os.path.join(output_path, f"blast_{cluster_type}_{id_}.tsv")
                # Write the dictionary to the file
                alignment_dict_to_file(write_dict, report_file_path, 'w', add_groups_ids)
        
        return add_groups_ids

    def process_cluster(class_,id_ , cluster, all_alleles, cds):
        """
        Process a single cluster.

        Parameters
        ----------
        class_ : str
            Class of the cluster.
        id_ : str or int
            ID of the cluster.
        cluster : str or int
            ID of the cluster.
        all_alleles : dict
            Dictionary of all alleles.
        cds : dict or str
            If single CDS then contain str if joined cluster then a dict.

        Returns
        -------
        id_ : str or int
            ID of the cluster.
        cluster : list
            List of the clusters.
        cluster_type : str
            Type of the cluster.
        is_cds : bool
            True if it's a CDS and not a loci, False otherwise.
        add_groups_ids : bool
            True if additional group IDs are present, False otherwise.

        """
        # Check the class and process accordingly
        if class_ == '1a':
            cluster_type = 'joined_cluster'
            cluster = cds[id_]
        else:
            id_ = cluster
            cluster = [cluster]
            cluster_type = 'retained'

        # Check if all_alleles exist
        if all_alleles:
            add_groups_ids = True
            is_cds = False
            cluster_alleles = []
            for entry in cluster:
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
            add_groups_ids = False
            is_cds = True

        return id_, cluster, cluster_type, is_cds, add_groups_ids

    def generate_write_dict(id_, cluster, is_cds, cds_matched_loci, representative_blast_results):
        """
        Generate the dictionary to be written to file.

        Parameters
        ----------
        id_ : str or int
            ID of the cluster.
        cluster : list
            List of the clusters.
        is_cds : bool
            Boolean indicating if it's a CDS.
        cds_matched_loci : dict
            Dictionary of CDS matched loci.
        representative_blast_results : dict
            Dictionary of representative blast results.

        Returns
        -------
        write_dict : dict
            Dictionary to be written to file.
        """
        # Check if it's a CDS and if it matches loci
        if is_cds and cds_matched_loci:
            queries = []
            if type(id_) == int:
                queries = cds_matched_loci[id_]
            else:
                for c in cluster:
                    queries += cds_matched_loci[c]
            # Generate the dictionary to be written
            write_dict = {query : {subject: {id_: entry for id_, entry in entries.items()}
                                for subject, entries in subjects.items() if subject in cluster}
                        for query, subjects in representative_blast_results.items()
                        if query in queries}
        else:
            # Generate the dictionary to be written
            write_dict = {query : {subject: {id_: entry for id_, entry in entries.items()}
                                for subject, entries in subjects.items()}
                        for query, subjects in representative_blast_results.items()
                        if query in cluster}
        return write_dict

    def process_classes(classes_outcome, representative_blast_results, output_path, add_groups_ids):
        """
        Process and write class results.

        Parameters
        ----------
        classes_outcome : list
            List of class outcomes.
        representative_blast_results : dict
            Dictionary of representative blast results.
        output_path : str
            Path to the output directory.
        add_groups_ids : bool
            Boolean indicating if additional group IDs are present.

        Returns
        -------
        No returns, writes files in output path.
        """
        # Loop over each class in the outcome
        for class_ in classes_outcome:
            # Generate the dictionary to be written
            write_dict = {query : {subject: {id_: entry for id_, entry in entries.items() if entry['class'] == class_}
                                for subject, entries in subjects.items()}
                        for query, subjects in representative_blast_results.items()}
            # Define the path of the report file
            report_file_path = os.path.join(output_path, f"blastn_group_{class_}.tsv")
            # Write the dictionary to the file
            alignment_dict_to_file(write_dict, report_file_path, 'w', add_groups_ids)

    # Create directories for output
    blast_by_cluster_output = os.path.join(output_path, 'blast_by_cluster')
    ff.create_directory(blast_by_cluster_output)
    blast_results_by_class_output = os.path.join(output_path, 'blast_results_by_class')
    ff.create_directory(blast_results_by_class_output)

    # Process and write cluster results
    add_groups_ids = process_clusters(cds_to_keep, representative_blast_results, all_alleles, cds_matched_loci, blast_by_cluster_output)

    # Process and write class results
    process_classes(classes_outcome, representative_blast_results, blast_results_by_class_output, add_groups_ids)

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
        if class_ in ['2b', '4b', '5b']:
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
        Translate possible new loci and write to master file.
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
            for cds in groups_paths_old:
                cds_name = f"retained_not_matched_by_blastn_{cds}"
                file_path = os.path.join(cds_outcome_results, cds_name)
                ff.copy_file(Retained_not_matched_by_blastn[cds], file_path)

        master_file_rep = None
        reps_trans_dict_cds = None
    else:
        # Initialize dictionaries to store paths
        groups_paths = {}
        groups_paths_reps = {}
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
    for column in ['Global_palign_all_min', 'Global_palign_all_max', 'Global_palign_pident_min', 'Global_palign_pident_max', 'Palign_local_min']:
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
    [cds_to_keep, important_relationships, drop_list] = process_classes(representative_blast_results,
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

    [_, _, _] = wrap_up_blast_results(cds_to_keep,
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
    [cds_to_keep, important_relationships, drop_list] = process_classes(representative_blast_results,
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
                                    representative_blast_results,
                                    classes_outcome,
                                    None,
                                    None,
                                    results_output)
    
    print("\nWrapping up BLAST results...")
    [groups_paths,
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