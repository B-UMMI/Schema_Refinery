import os
import concurrent.futures
from itertools import repeat

try:
    from utils import (file_functions as ff,
                       sequence_functions as sf,
                       clustering_functions as cf,
                       blast_functions as bf,
                       alignments_functions as af,
                       iterable_functions as itf,
                       linux_functions as lf,
                       graphical_functions as gf,
                       kmers_functions as kf,
                       pandas_functions as pf)
except ModuleNotFoundError:
    from SchemaRefinery.utils import (file_functions as ff,
                                      sequence_functions as sf,
                                      clustering_functions as cf,
                                      blast_functions as bf,
                                      alignments_functions as af,
                                      iterable_functions as itf,
                                      linux_functions as lf,
                                      graphical_functions as gf,
                                      kmers_functions as kf,
                                      pandas_functions as pf)

def alignment_dict_to_file(blast_results_dict, file_path, write_type, add_group_column = False):
    """
    Writes alignment data to a file from a nested dictionary structure.

    This function takes a nested dictionary containing alignment data, a file path, and a write type (write or append)
    as input and writes the alignment data to the specified file. It supports the option to add an additional column
    for CDS group information in the output file's header.

    Parameters
    ----------
    blast_results_dict : dict
        A nested dictionary where the first level keys are query IDs, the second level keys are subject IDs, and the
        third level keys are specific alignment IDs, mapping to dictionaries containing alignment data.
    file_path : str
        The path to the file where the alignment data should be written or appended.
    write_type : str
        Specifies whether to create a new file ('w') and write the data or to append ('a') the data to an existing file.
    add_group_column : bool, optional
        Indicates whether to add a 'CDS_group' column to the header of the output file. Defaults to False.

    Returns
    -------
    None
        Writes the alignment data to the specified file based on the provided dictionary structure and parameters.

    Notes
    -----
    - The function constructs a header for the output file based on the alignment data structure and the `add_group_column`
    parameter.
    - It iterates over the nested dictionary structure to write each piece of alignment data to the file, formatting
    each row as tab-separated values.
    - This function is useful for exporting BLAST alignment results to a file for further analysis or reporting.
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
              'Frequency_in_genomes_query\t',
              'Frequency_in_genomes_subject\t',
              'Global_palign_all_min\t',
              'Global_palign_all_max\t',
              'Global_palign_pident_min\t',
              'Global_palign_pident_max\t',
              'Palign_local_min\t',
              'Class\n']

    if add_group_column:
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
                         frequency_in_genomes, loci_ids, add_groups_ids = None):
    """
    Enhances BLAST results with additional metrics and frequencies.

    This function enriches the given BLAST results dictionary with several key metrics and frequencies
    to provide a more comprehensive analysis of the BLAST hits. It adds metrics such as Blast
    Score Ratio (BSR), k-mer similarities and coverage, and frequencies of query and subject CDS in
    schema genomes. It also includes global and local pairwise alignment scores, both in terms of
    coverage and percentage identity, with the ability to focus on specific percentage identity thresholds.

                                  
    Parameters
    ----------                                  
    representative_blast_results : dict
        A dictionary containing BLAST results. Each entry is expected to represent a unique BLAST hit with
        various metrics.
    reps_kmers_sim : dict
        A dictionary mapping pairs of CDS to their k-mer similarity scores.
    bsr_values : dict
        A dictionary mapping pairs of CDS to their Blast Score Ratio (BSR) values. Can be None if BSR values
        are not available.
    representative_blast_results_coords_all : dict
        A dictionary containing the coordinates for all BLAST entries, used for calculating global alignment metrics.
    representative_blast_results_coords_pident : dict
        A dictionary containing the coordinates for BLAST entries above a certain percentage identity threshold,
        used for calculating specific global alignment metrics.
    frequency_in_genomes : dict
        A dictionary summarizing the frequency of each representative cluster within the genomes of the schema,
        enhancing the context of BLAST results.
    loci_ids : bool
        Indicates whether the IDs of loci representatives are included in the `frequency_in_genomes`. If true,
        IDs follow the format `loci1_x`.
    add_groups_ids : Dict, optional
        A dictionary mapping group IDs to their member CDS. This is used to add group information to the BLAST
        results for enhanced analysis.

    Returns
    -------
    No returns, modifies the representative_blast_results dict inside the main
    function.

    Notes
    -----
    - The function is designed to work with detailed BLAST results and requires several pre-computed metrics
    and frequencies as input.
    - It is crucial for enhancing the analysis of BLAST results, especially in comparative genomics and
    schema development projects.
    """
    def get_kmer_values(reps_kmers_sim, query, subject):
        """
        Retrieves k-mer similarity and coverage values for a specified query and subject pair.

        This function looks up the k-mer similarity and coverage values between a query and a
        subject sequence from a precomputed dictionary. It is designed to facilitate the analysis
        of genomic sequences by providing key metrics that reflect the degree of similarity and the
        extent of coverage between pairs of sequences. These metrics are crucial for understanding
        the genetic relationships and variations among sequences.

        Parameters
        ----------
        reps_kmers_sim : dict
            A dictionary where keys are query sequence IDs and values are dictionaries with subject
            sequence IDs as keys. Each inner dictionary's values are tuples containing the k-mer
            similarity and coverage values.
        query : str
            The identifier for the query sequence. This is used to look up the corresponding dictionary
            of subjects within `reps_kmers_sim`.
        subject : str
            The identifier for the subject sequence. This is used to retrieve the similarity and coverage
            values from the dictionary associated with the `query`.

        Returns
        -------
        sim : float or str
            The k-mer similarity value between the query and subject sequences. Returns 0 if no value is
            found, or '-' if the `reps_kmers_sim` dictionary is empty or not provided.
        cov : float or str
            The coverage value indicating the extent to which the k-mers of the query sequence are present
            in the subject sequence. Returns 0 if no value is found, or '-' if the `reps_kmers_sim` dictionary
            is empty or not provided.

        Notes
        -----
        - The function is robust to missing data, providing default values when specific similarity or coverage values
        are unavailable.
        - It is a utility function primarily used in genomic analysis workflows, particularly in the context of
        comparing sequences for similarity and coverage using k-mer based metrics.
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
        Fetches the BLAST Score Ratio (BSR) for a specified pair of sequences.

        This function extracts the BSR value for a given pair of sequences identified by their respective
        query and subject IDs from a pre-populated dictionary. The BSR is a normalized metric used to
        compare the BLAST scores of different alignments, providing insight into the relative similarity
        between sequences. If the BSR exceeds 1.0, it is rounded to the nearest whole number to maintain
        consistency in reporting.

        Parameters
        ----------
        bsr_values : dict
            A dictionary where keys are query sequence IDs and values are dictionaries with subject sequence
            IDs as keys. Each inner dictionary's values are the BSR values.
        query : str
            The identifier for the query sequence. This is used to select the appropriate dictionary of subject
            sequences within `bsr_values`.
        subject : str
            The identifier for the subject sequence. This is used to retrieve the BSR value from the dictionary
            associated with the `query`.

        Returns
        -------
        bsr : float
            The BSR value between the query and subject sequences. Returns 0 if no BSR value is found for the
            given pair. If the BSR value is greater than 1, it is rounded to the nearest whole number. The BSR
            value is rounded by 4 decimal places for consistency.

        Notes
        -----
        - The BSR value is a crucial metric in bioinformatics for assessing the quality of sequence alignments, with values typically ranging from 0 to 1. Values greater than 1 are considered anomalies and are adjusted accordingly.
        - This function is essential for workflows involving comparative genomics or sequence alignment analysis, where BSR values provide a standardized measure of sequence similarity.
        """
        bsr = bsr_values[query].get(subject, 0)
        if bsr > 1.0:
            bsr = float(round(bsr))
        return round(bsr, 4)

    def calculate_total_length(representative_blast_results_coords, query, subject):
        """
        Calculates the total aligned length for each reference sequence in a given query-subject pair.

        This function computes the total length of aligned sequences for each reference sequence associated
        with a specific query-subject pair. It processes the alignment intervals for each reference sequence,
        merges overlapping intervals to avoid double counting, and sums up the lengths of these intervals to
        determine the total aligned length.

        Parameters
        ----------
        representative_blast_results_coords : dict
            A nested dictionary where the first level keys are query sequence IDs, the second level keys are
            subject sequence IDs, and the values are dictionaries mapping reference sequence IDs to lists of
            alignment intervals.
        query : str
            The identifier for the query sequence. This is used to select the appropriate dictionary of subject
            sequences within `representative_blast_results_coords`.
        subject : str
            The identifier for the subject sequence. This is used to retrieve the dictionary of reference sequences
            and their alignment intervals.

        Returns
        -------
        total_length : dict
            A dictionary where keys are reference sequence IDs and values are the total aligned length for that
            reference sequence. The total aligned length is calculated by merging overlapping intervals and
            summing the lengths of the resulting intervals.

        Notes
        -----
        - The function assumes that alignment intervals are provided as tuples or lists of two elements, where the
        first element is the start position and the second element is the end position of the interval.
        - Overlapping intervals for each reference sequence are merged to ensure that the total aligned length is
        accurately calculated without double counting overlapping regions.
        - This function is particularly useful in genomic analyses where understanding the extent of alignment
        coverage is important for interpreting BLAST results.
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
        Calculates the minimum and maximum global pairwise alignment percentages.

        This function computes the global pairwise alignment (palign) percentages for a query and its subject
        based on their total aligned lengths and original lengths. It calculates both the minimum and maximum
        palign values to provide a range of alignment coverage, which can be useful for assessing the quality
        and extent of the alignment between the query and subject sequences.

        Parameters
        ----------
        total_length : dict
            A dictionary where keys are 'query' and 'subject', and values are the total aligned lengths for the
            query and subject sequences, respectively.
        result : dict
            A dictionary containing the lengths of the query and subject sequences under the keys 'query_length'
            and 'subject_length'.

        Returns
        -------
        global_palign_min : float
            The minimum global pairwise alignment percentage, calculated as the smaller of the two ratios: total
            aligned length of the query to its original length, and total aligned length of the subject to its
            original length. The value is rounded to 4 decimal places.
        global_palign_max : float
            The maximum global pairwise alignment percentage, calculated as the larger of the two ratios: total
            aligned length of the query to its original length, and total aligned length of the subject to its
            original length. The value is rounded to 4 decimal places.

        Notes
        -----
        - The global pairwise alignment percentage is a measure of how much of the original sequences
        (query and subject) are covered by the alignment. It provides insight into the completeness of the alignment.
        - This function is particularly useful in bioinformatics for evaluating the quality of sequence alignments,
        where higher coverage percentages might indicate more reliable alignments.
        """
        global_palign_min = min(total_length['query'] / result['query_length'],
                                total_length['subject'] / result['subject_length'])
        global_palign_max = max(total_length['query'] / result['query_length'],
                                total_length['subject'] / result['subject_length'])
        return round(global_palign_min, 4), round(global_palign_max, 4)

    def calculate_local_palign(result):
        """
        Calculates the minimum local pairwise alignment percentage.

        This function computes the local pairwise alignment (palign) percentage for a given BLAST search result
        by comparing the aligned lengths of the query and subject sequences to their total lengths. It calculates
        the alignment percentage for both the query and subject, and returns the minimum of these two percentages.
        This metric is useful for assessing the extent of alignment within the local regions of interest in both 
        sequences.

        Parameters
        ----------
        result : dict
            A dictionary containing the result of a BLAST search, including the start and end positions of the alignment
            on both the query and subject sequences, as well as their total lengths.

        Returns
        -------
        local_palign_min : float
            The minimum local pairwise alignment percentage, calculated as the smaller of the two ratios: aligned
            length of the query to its total length, and aligned length of the subject to its total length. The value
            is rounded to 4 decimal places.

        Notes
        -----
        - The local pairwise alignment percentage provides insight into the local similarity between the query and
        subject sequences, focusing on the aligned regions.
        - This function is particularly useful in sequence alignment analyses, where understanding the coverage and
        similarity of local alignments is important.
        """
        local_palign_min = min((result['query_end'] - result['query_start'] + 1) / result['query_length'],
                            (result['subject_end'] - result['subject_start'] + 1) / result['subject_length'])
        return round(local_palign_min, 4)

    def update_results(representative_blast_results, query, subject, entry_id, bsr, sim, cov, frequency_in_genomes,
                       global_palign_all_min, global_palign_all_max, global_palign_pident_min, global_palign_pident_max,
                       local_palign_min, loci_ids, add_groups_ids):
        """
        Updates the BLAST results for a specific query and subject pair with new data.

        This function modifies the existing BLAST results dictionary by updating the entries for a given query and
        subject pair with new information. It handles the addition of various metrics such as BSR, similarity,
        coverage, frequency in genomes, global and local pairwise alignment percentages, and group IDs. The function
        also supports the modification of query and subject IDs based on loci information.

        Parameters
        ----------
        representative_blast_results : dict
            A dictionary containing BLAST results where keys are query IDs, values are dictionaries with subject
            IDs as keys, and each subject dictionary contains dictionaries of entry IDs with their respective data.
        query : str
            The query sequence ID.
        subject : str
            The subject sequence ID.
        entry_id : str
            The ID of the entry to update within the BLAST results.
        bsr, sim, cov : float
            The BSR, similarity, and coverage values to update.
        frequency_in_genomes : dict
            A dictionary containing the frequency of the query and subject in genomes.
        global_palign_all_min, global_palign_all_max, global_palign_pident_min, global_palign_pident_max : float
            The minimum and maximum global pairwise alignment percentages, including those based on Pident threshold.
        local_palign_min : float
            The minimum local pairwise alignment percentage.
        loci_ids : list
            A list indicating whether to modify the query and/or subject IDs based on loci information.
        add_groups_ids : dict
            A dictionary containing group IDs to be added to the results, where keys are subject IDs and values
            are the group members.

        Returns
        -------
        None
            This function does not return any value but modifies the `representative_blast_results` dictionary in place.

        Notes
        -----
        - The function assumes the presence of a utility module `itf` with functions for ID manipulation and
        identification within dictionaries.
        - It is designed to be flexible, allowing for the update of specific metrics as needed without requiring
        a complete overhaul of the entry data.
        """
        if loci_ids[0]:
            query_before = query
            query = itf.remove_by_regex(query, '_(\d+)')
        if loci_ids[1]:
            subject_before = subject
            subject = itf.remove_by_regex(subject, '_(\d+)')
        update_dict = {
            'bsr': bsr,
            'kmers_sim': sim,
            'kmers_cov': cov,
            'frequency_in_genomes_query_cds': frequency_in_genomes[query],
            'frequency_in_genomes_subject_cds': frequency_in_genomes[subject],
            'global_palign_all_min' : global_palign_all_min,
            'global_palign_all_max': global_palign_all_max,
            'global_palign_pident_min': global_palign_pident_min,
            'global_palign_pident_max': global_palign_pident_max,
            'local_palign_min': local_palign_min
        }
        if loci_ids[0]:
            query = query_before
        if loci_ids[1]:
            subject = subject_before
        representative_blast_results[query][subject][entry_id].update(update_dict)

        if add_groups_ids:
            id_ = itf.identify_string_in_dict(subject, add_groups_ids)
            if not id_:
                id_ = subject
            update_dict = {'cds_group': id_}
            representative_blast_results[query][subject][entry_id].update(update_dict)

    def remove_results(representative_blast_results, query, subject, entry_id):
        """
        Removes a specific entry from the BLAST results for a given query and subject pair.

        This function is designed to modify an existing dictionary of BLAST results by removing a specified
        entry identified by its entry ID for a particular query and subject pair. It directly alters the
        `representative_blast_results` dictionary, removing the entry corresponding to the provided `entry_id`
        within the nested structure of query and subject IDs.

        Parameters
        ----------
        representative_blast_results : dict
            A dictionary containing BLAST results, structured with query IDs as keys, each mapping to a dictionary
            of subject IDs, which in turn map to dictionaries of entry IDs and their associated data.
        query : str
            The identifier for the query sequence, used to locate the correct subset of results within
            `representative_blast_results`.
        subject : str
            The identifier for the subject sequence, used in conjunction with `query` to further narrow down the
            specific subset of results.
        entry_id : str
            The identifier of the specific entry to be removed from the results.

        Returns
        -------
        None
            This function does not return any value. It modifies the `representative_blast_results` dictionary in
            place, removing the specified entry.

        Notes
        -----
        - This function is useful for cleaning up BLAST results, allowing for the removal of specific entries that are
        no longer needed or relevant.
        - It operates directly on the provided dictionary, requiring careful handling to avoid unintended modifications.
        """
        del representative_blast_results[query][subject][entry_id]

    def clean_up_results(representative_blast_results, query, subject):
        """
        Cleans up BLAST results for a specific query and subject by removing empty entries.

        This function is designed to modify an existing dictionary of BLAST results by checking for and
        removing any entries that are empty for a given query and subject pair. It aims to streamline the BLAST
        results by ensuring that only entries with data are retained. The function operates directly on the
        `representative_blast_results` dictionary, removing entries without returning any value.

        Parameters
        ----------
        representative_blast_results : dict
            A dictionary containing BLAST results, where keys are query IDs, and values are dictionaries with
            subject IDs as keys, each mapping to their respective result entries.
        query : str
            The identifier for the query sequence, used to locate the correct subset of results within
            `representative_blast_results`.
        subject : str
            The identifier for the subject sequence, used in conjunction with `query` to further narrow down the
            specific subset of results to be cleaned.
        
        Returns
        -------
        None
            This function does not return any value. It modifies the `representative_blast_results` dictionary in
            place, removing the specified entry.

        Notes
        -----
        - The function checks for and removes entries that are empty for the specified query and subject. If the
        subject entry under a query is empty, it is removed. If this results in the query entry becoming empty,
        it is also removed.
        - This cleanup process is essential for maintaining the integrity and usability of BLAST results, especially
        in large-scale genomic analyses where empty entries can clutter the dataset.
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

            total_length = calculate_total_length(representative_blast_results_coords_all, query, subject)
            global_palign_all_min, global_palign_all_max = calculate_global_palign(total_length, blastn_results[1])

            total_length = calculate_total_length(representative_blast_results_coords_pident, query, subject)
            global_palign_pident_min, global_palign_pident_max = calculate_global_palign(total_length, blastn_results[1])
            # Iterate over the blastn_results dictionary
            for entry_id, result in list(blastn_results.items()):

                local_palign_min = calculate_local_palign(result)
                # Remove entries with negative local palign values meaning that they are inverse alignments.
                if local_palign_min >= 0:
                    update_results(representative_blast_results, query, subject, entry_id, bsr, sim, cov, frequency_in_genomes, global_palign_all_min, global_palign_all_max, global_palign_pident_min, global_palign_pident_max, local_palign_min, loci_ids, add_groups_ids)
                else:
                    remove_results(representative_blast_results, query, subject, entry_id)

            clean_up_results(representative_blast_results, query, subject)

def separate_blastn_results_into_classes(representative_blast_results, constants):
    """
    Separates BLAST results into predefined classes based on specific criteria.

    This function iterates through BLAST results and classifies each result into a specific class
    based on criteria such as global alignment percentage, bit score ratio (bsr), and frequency ratios
    between query and subject CDS in genomes. The classification is done by updating the results
    dictionary with a new key-value pair indicating the class of each BLAST result.

    Parameters
    ----------
    representative_blast_results : dict
        A nested dictionary where the first level keys are query sequence IDs, the second level keys
        are subject sequence IDs, and the third level keys are unique identifiers for each BLAST
        result. Each BLAST result is a dictionary containing keys such as 'frequency_in_genomes_query_cds',
        'frequency_in_genomes_subject_cds', 'global_palign_all_min', 'bsr', and 'pident'.
    constants : tuple or list
        A collection of constants used in the classification criteria. Specifically, `constants[1]`
        is used as a threshold for the percentage identity (pident) in one of the classification conditions.

    Returns
    -------
    classes_outcome : tuple
        A tuple of class identifiers indicating the order of priority for the classes.

    Notes
    -----
    - The function modifies `representative_blast_results` in place by adding a 'class' key to each
      BLASTN result dictionary.
    - The classification logic is based on a combination of alignment quality metrics and frequency
      ratios, with specific thresholds and conditions determining the class assignment.
    - The function assumes that `constants` provides necessary thresholds for classification and
      that its elements are accessed by index.
    """
    def add_class_to_dict(class_name):
        """
        Adds a class identifier to a BLAST result within the representative_blast_results dictionary.

        This helper function is used to update the BLASTNresult dictionaries with a 'class' key,
        assigning the specified class identifier based on the classification logic in the outer function.

        Parameters
        ----------
        class_name : str
            The class identifier to be added to the BLASTN result. This should be one of the values
            from the classes_outcome tuple defined in the outer function.

        Notes
        -----
        - This function directly modifies the `representative_blast_results` dictionary from the outer
          scope, specifically adding or updating the 'class' key for a BLASTN result.
        - It is designed to be used only within the `separate_blastn_results_into_classes` function.
        """
        representative_blast_results[query][id_subject][id_].update({'class': class_name})

    # Define classes based on priority
    classes_outcome = ('1a', '1b', '2a', '3a', '2b', '1c', '3b', '4a', '4b', '4c','5')

    # Loop through the representative BLAST results
    for query, rep_blast_result in representative_blast_results.items():
        for id_subject, matches in rep_blast_result.items():
            for id_, blastn_entry in matches.items():
                # Calculate the frequency ratio
                query_freq = blastn_entry['frequency_in_genomes_query_cds']
                subject_freq = blastn_entry['frequency_in_genomes_subject_cds']
                # If one of the frequencies is 0, set the ratio to 0.1 if the other frequency is 10 times greater
                if query_freq == 0 or subject_freq == 0:
                    freq_ratio = 0.1 if query_freq > 10 or subject_freq > 10 else 1
                else:
                    freq_ratio = min(query_freq/subject_freq,
                                    subject_freq/query_freq)
                
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

def sort_blast_results_by_classes(representative_blast_results, classes_outcome):
    """
    Sorts BLAST results by classes based on the alignment score.

    This function organizes BLAST results into a sorted structure according to predefined classes.
    It ensures that for each query, the results are grouped by the class of the alignment, prioritizing
    the classes as specified in the `classes_outcome` list.

    Parameters
    ----------
    representative_blast_results : dict
        A dictionary where each key is a query identifier and each value is another dictionary.
        The inner dictionary's keys are subject identifiers, and values are lists containing
        details of the match, where the second element is a dictionary with the key 'class'
        indicating the class of the alignment.
    classes_outcome : tuple
        A list of possible classes outcomes to sort the BLAST results into. The order in this list
        determines the priority of the classes when organizing the results.

    Returns
    -------
    sorted_blast_dict : dict
        A dictionary structured similarly to `representative_blast_results`, but sorted such that
        all results for a given query are grouped by their class as determined by the highest
        scoring alignment.

    Notes
    -----
    - The function assumes that each match list in the values of `representative_blast_results`
      contains at least one element, which is a dictionary with a 'class' key.
    - It creates a temporary dictionary to first group results by class, then consolidates these
      into the final sorted dictionary to be returned.
    """
    sorted_blast_dict = {}
    temp_dict = {k: {} for k in classes_outcome}
    
    for query, rep_blast_result in representative_blast_results.items():
        for id_subject, matches in rep_blast_result.items():
            # Class of the alignment with biggest score.
            class_ = matches[1]['class']
            if not temp_dict[class_].get(query):
                temp_dict[class_][query] = {}
            temp_dict[class_][query][id_subject] = matches
    
    for class_, sorted_blast_reps in temp_dict.items():
        for query, rep_blast_result in sorted_blast_reps.items():
            if not sorted_blast_dict.get(query):
                sorted_blast_dict[query] = {}
            for id_subject, matches in rep_blast_result.items():
                    sorted_blast_dict[query][id_subject] = matches
    
    return sorted_blast_dict

def process_classes(representative_blast_results, classes_outcome, all_alleles = None):
    """
    Processes BLAST results to determine class-based relationships and counts.

    This function iterates through representative BLAST results to establish relationships
    between different coding sequences (CDS) and to count occurrences by class. It handles
    allele replacements, prioritizes classes based on a predefined order, and identifies
    important relationships between sequences.

    Parameters
    ----------
    representative_blast_results : dict
        A nested dictionary where the first key is the query sequence ID, the second key is
        the subject sequence ID, and the value is another dictionary containing match details
        including the class of the match.
    classes_outcome : tuple
        A list of class identifiers ordered by priority. This order determines which classes are
        considered more significant when multiple matches for the same pair of sequences are found.
    all_alleles : dict, optional
        A dictionary mapping sequence IDs to their corresponding allele names. If provided, it is
        used to replace allele IDs with loci/CDS names in the processing.

    Returns
    -------
    processed_results : dict
        A dictionary containing processed results with keys formatted as "query|subject" and values being tuples
        containing information about the processed sequences, their class, relationships, and additional details.
    count_results_by_class : dict
        A dictionary containing counts of results by class, with keys formatted as "query|subject" and values being
        dictionaries with class identifiers as keys and counts as values.
    reps_and_alleles_ids : dict
        A dictionary mapping pairs of query and subject sequences to their unique loci/CDS IDs and alleles IDs.

    Notes
    -----
    - The function dynamically adjusts based on the presence of `all_alleles`, affecting how sequence IDs
    are replaced and processed.
    - It employs a complex logic to handle different scenarios based on class types and the presence or absence of
    alleles in the processed results, including handling allele replacements and determining the importance of
    relationships.
    """
    # Initialize variables
    count_results_by_class = {}
    reps_and_alleles_ids = {}
    processed_results = {}
    drop_mark = []
    # Process the CDS to find what CDS to retain while also adding the relationships between different CDS
    for query, rep_blast_result in representative_blast_results.items():
        for id_subject, matches in rep_blast_result.items():
            class_ = matches[1]['class'] 
            ids_for_relationship = [query, id_subject]
            new_query = query
            new_id_subject = id_subject

            strings = [str(query), str(id_subject), class_]
            if all_alleles:
                replaced_query = itf.identify_string_in_dict(query, all_alleles)
                if replaced_query:
                    new_query = replaced_query
                    strings[0] = new_query if isinstance(new_query, str) else f"{query}({new_query})"
                replaced_id_subject = itf.identify_string_in_dict(id_subject, all_alleles)
                if replaced_id_subject:
                    new_id_subject = replaced_id_subject
                    strings[1] = new_id_subject if isinstance(new_id_subject, str) else f"{id_subject}({new_id_subject})"

                current_allele_class_index = classes_outcome.index(class_)
                # Check if the current loci were already processed
                if not processed_results.get(f"{new_query}|{new_id_subject}"):
                    run_next_step = True
                # If those loci/CDS were already processed, check if the current class is better than the previous one
                elif current_allele_class_index < classes_outcome.index(processed_results[f"{new_query}|{new_id_subject}"][0]):
                    run_next_step = True
                # If not then skip the current alleles
                else:
                    run_next_step = False
            else:
                run_next_step = True

            count_results_by_class.setdefault(f"{new_query}|{new_id_subject}", {})
            if not count_results_by_class[f"{new_query}|{new_id_subject}"].get(class_):
                count_results_by_class[f"{new_query}|{new_id_subject}"].setdefault(class_, 1)
            else:
                count_results_by_class[f"{new_query}|{new_id_subject}"][class_] += 1
            # Get unique loci/CDS for each query and subject rep and allele.
            reps_and_alleles_ids.setdefault(f"{new_query}|{new_id_subject}", [set(), set()])
            if ids_for_relationship[0] not in reps_and_alleles_ids[f"{new_query}|{new_id_subject}"][0]:
                reps_and_alleles_ids[f"{new_query}|{new_id_subject}"][0].add(ids_for_relationship[0])
            if ids_for_relationship[1] not in reps_and_alleles_ids[f"{new_query}|{new_id_subject}"][1]:
                reps_and_alleles_ids[f"{new_query}|{new_id_subject}"][1].add(ids_for_relationship[1])
    
            if run_next_step:
                # Set all None to run newly for this query/subject combination
                processed_results[f"{new_query}|{new_id_subject}"] = (None,
                                                        None,
                                                        None,
                                                        None,
                                                        None,
                                                        None)

                if class_ in ['1b', '2a', '3a']:
                    blastn_entry = matches[list(matches.keys())[0]]
                    # Determine if the frequency of the query is greater than the subject.
                    is_frequency_greater = blastn_entry['frequency_in_genomes_query_cds'] >= blastn_entry['frequency_in_genomes_subject_cds']
                    # Determine if the query or subject should be dropped.
                    query_or_subject = new_id_subject if is_frequency_greater else new_query
                else:
                    query_or_subject = []

                # For the related_matches.tsv file.
                if class_ not in ['4c','5']:
                    # Add asterisk to the query or subject that was dropped.
                    if class_ in ['1b', '2a', '3a']:
                        blastn_entry = matches[list(matches.keys())[0]]
                        # Determine if the frequency of the query is greater than the subject.
                        is_frequency_greater = blastn_entry['frequency_in_genomes_query_cds'] >= blastn_entry['frequency_in_genomes_subject_cds']
                        # Determine if the query or subject should be dropped.
                        dropped = new_id_subject if is_frequency_greater else new_query
                        if new_query == dropped:
                            drop_mark.append(new_query)
                            strings[0] += '*' 
                        else:
                            drop_mark.append(new_id_subject)
                            strings[1] += '*'

                processed_results[f"{new_query}|{new_id_subject}"] = (class_,
                                                    ids_for_relationship,
                                                    query_or_subject,
                                                    (new_query, new_id_subject),
                                                    strings)

    return processed_results, count_results_by_class, reps_and_alleles_ids, drop_mark

def extract_results(processed_results, count_results_by_class, frequency_in_genomes,
                    cds_to_keep, drop_set, run_normal, classes_outcome):
    """
    Extracts and organizes results from process_classes.

    Parameters
    ----------
    processed_results : dict
        The processed results data.
    count_results_by_class : dict
        A dictionary with counts of results by class.
    frequency_in_genomes : dict
        A dictionary containing the frequency of the query and subject in genomes.
    classes_outcome : list
        A list of class outcomes.

    Returns
    -------
    all_relationships : dict
        All relationships between loci and CDS.
    related_clusters : dict
        Dict that groups CDS/loci by ID and that contains strings to write in output file.
    
    Notes
    -----
    - The function iterates over `processed_results` to organize and cluster related CDS/loci based
    on their classification outcomes and the presence in specific clusters.
    - It uses helper functions like `cf.cluster_by_ids` for clustering and `itf.identify_string_in_dict`
    for identifying if a query or subject ID is present in the clusters.
    """
    def cluster_data(run_normal, processed_results):
        if run_normal:
            key_extractor = lambda v: v[3]
            condition = lambda v: v[0] not in ['4c', '5']
        else:
            key_extractor = lambda v: [itf.remove_by_regex(v[4][0], regex_pattern), itf.remove_by_regex(v[4][1], regex_pattern)]
            condition = lambda v: v[0] not in ['4c', '5']

        return {i: cluster for i, cluster in enumerate(cf.cluster_by_ids([key_extractor(v) for v in processed_results.values() if condition(v)]), 1)}

    def choice_data(run_normal, processed_results, to_cluster_list):
        if run_normal:
            key_extractor = lambda v: v[3]
            additional_condition = lambda v: '*' in v[4][0] or '*' in v[4][1]
        else:
            key_extractor = lambda v: [itf.remove_by_regex(v[4][0], regex_pattern), itf.remove_by_regex(v[4][1], regex_pattern)]
            additional_condition = lambda v: False  # No additional condition in non-normal run
    
        return {i: cluster for i, cluster in enumerate(cf.cluster_by_ids([key_extractor(v) for v in processed_results.values() if v[0] in ['1c', '2b', '3b', '4b'] or ((itf.identify_string_in_dict(v[3][0], to_cluster_list) and additional_condition(v)) or (itf.identify_string_in_dict(v[3][1], to_cluster_list) and additional_condition(v)))]), 1)}
    
    def process_id(id_, run_normal, regex_pattern, to_cluster_list, cds_to_keep):
        processed_id = id_ if run_normal else itf.remove_by_regex(id_, regex_pattern)
        present = itf.identify_string_in_dict(processed_id, to_cluster_list)
        cds_joined_id = None if run_normal else id_
        joined_id = itf.identify_string_in_dict(processed_id if run_normal else cds_joined_id, cds_to_keep['1a'])
        return processed_id, present, cds_joined_id, joined_id

    def check_in_recommendations(id_, joined_id, recommendations, key, categories):
        return any((joined_id or id_) in itf.flatten_list([v for k, v in recommendations[key].items() if cat in k]) for cat in categories)

    def format_id(id_, cds_joined_id, run_normal):
        return id_ if run_normal or not isinstance(cds_joined_id, int) else f"{id_}({cds_joined_id})"

    def add_to_recommendations(category, id_to_write, joined_id=None):
        if joined_id is not None:  # For joined or choice categories
            recommendations[key].setdefault(f'{category}_{joined_id}', set()).add(id_to_write)
        else:  # For keep or drop categories
            recommendations[key].setdefault(category, set()).add(id_to_write)

    all_relationships = {class_: [] for class_ in classes_outcome}
    related_clusters = {}
    recommendations = {}
    processed_cases = []
    regex_pattern = r'\([^()]*\)|\*'
    # Normal run, where IDs are only loci or CDS original IDs.
    to_cluster_list = cluster_data(run_normal, processed_results)
    choice = choice_data(run_normal, processed_results, to_cluster_list)

    related_clusters = {}
    for results in processed_results.values():
        if results[0] in ['4c','5']:
            continue

        query_id, query_present, query_cds_joined_id, joined_query_id = process_id(results[3][0] if run_normal else results[4][0], run_normal, regex_pattern, to_cluster_list, cds_to_keep)
        subject_id, subject_present, subject_cds_joined_id, joined_subject_id = process_id(results[3][1] if run_normal else results[4][1], run_normal, regex_pattern, to_cluster_list, cds_to_keep)

        key = query_present if query_present else subject_present

        related_clusters.setdefault(key, []).append(results[4] 
                                                    + [f"{count_results_by_class[f'{results[3][0]}|{results[3][1]}'][results[0]]}/{sum(count_results_by_class[f'{results[3][0]}|{results[3][1]}'].values())}"]
                                                    + [str(frequency_in_genomes[results[3][0]])]
                                                    + [str(frequency_in_genomes[results[3][1]])])

        recommendations.setdefault(key, {})
        if_same_joined = (joined_query_id == joined_subject_id) if joined_query_id and joined_subject_id else False
        if_joined_query = check_in_recommendations(query_id, joined_query_id, recommendations, key, ['Joined'])
        if_joined_subject = check_in_recommendations(subject_id, joined_subject_id, recommendations, key, ['Joined'])
        if_query_in_choice = check_in_recommendations(query_id, joined_query_id, recommendations, key, ['Choice'])
        if_subject_in_choice = check_in_recommendations(subject_id, joined_subject_id, recommendations, key, ['Choice'])
        if_query_in_keep = check_in_recommendations(query_id, joined_query_id, recommendations, key, ['Keep'])
        if_subject_in_keep = check_in_recommendations(subject_id, joined_subject_id, recommendations, key, ['Keep'])
        if_query_dropped = (joined_query_id or query_id) in drop_set
        if_subject_dropped = (joined_subject_id or subject_id) in drop_set

        choice_query_id = itf.identify_string_in_dict(query_id, choice)
        choice_subject_id = itf.identify_string_in_dict(subject_id, choice)

        # What IDs to addto the Keep, Drop and Choice.
        query_to_write = format_id(joined_query_id or query_id, query_cds_joined_id, run_normal)
        subject_to_write = format_id(joined_subject_id or subject_id, subject_cds_joined_id, run_normal)
        #TODO teste this code
        query_to_write = query_to_write if not isinstance(joined_query_id, int) else f"Joined_{query_to_write}"
        subject_to_write = subject_to_write if not isinstance(joined_subject_id, int) else f"Joined_{subject_to_write}"

        joined_query_to_write = format_id(query_id, query_cds_joined_id, run_normal)
        joined_subject_to_write = format_id(subject_id, subject_cds_joined_id, run_normal)

        # Check if the pair was not processed yet
        if [query_id, subject_id] not in processed_cases:
            processed_cases.append([subject_id, query_id])  # Add the inverse pair to the processed cases

            if results[0] == '1a':
                if isinstance(joined_query_id, int):
                    add_to_recommendations('Joined', joined_query_to_write, joined_query_id)
                if isinstance(joined_subject_id, int):
                    add_to_recommendations('Joined', joined_subject_to_write, joined_subject_id)

            elif results[0] in ['1c', '2b', '3b', '4b']:
                if if_query_in_keep:
                    recommendations[key]['Keep'].remove(query_to_write)
                    if len(recommendations[key]['Keep']) == 0:
                        recommendations[key].pop('Keep')
                if if_subject_in_keep:
                    recommendations[key]['Keep'].remove(subject_to_write)
                    if len(recommendations[key]['Keep']) == 0:
                        recommendations[key].pop('Keep')
                if not if_query_dropped and not if_subject_dropped and not if_same_joined:
                    add_to_recommendations('Choice', query_to_write, choice_query_id)
                    add_to_recommendations('Choice', subject_to_write, choice_subject_id)

            elif results[0] in ['1b', '2a', '3a', '4a']:
                if (joined_query_id and '*' in results[4][0]) or (joined_subject_id and '*' in results[4][1]):
                    add_to_recommendations('Choice', query_to_write, choice_query_id)
                    add_to_recommendations('Choice', subject_to_write, choice_subject_id)
                if query_id in drop_set:
                    if not if_joined_subject and not if_subject_in_choice:
                        add_to_recommendations('Keep', subject_to_write)
                    if not if_joined_query and not if_query_in_choice:
                        add_to_recommendations('Drop', query_to_write)
                elif subject_id in drop_set:
                    if not if_joined_query and not if_query_in_choice:
                        add_to_recommendations('Keep', query_to_write)
                    if not if_joined_subject and not if_subject_in_choice:
                        add_to_recommendations('Drop', subject_to_write)

    for k, v in processed_results.items():
        all_relationships.setdefault(v[0], []).append(v[1])

    sort_order = ['Joined', 'Choice', 'Keep', 'Drop']
    recommendations = {k: {l[0]: l[1] for l in sorted(v.items(), key=lambda x: sort_order.index(x[0].split('_')[0]))} for k, v in recommendations.items()}
    
    return all_relationships, related_clusters, recommendations

def write_blast_summary_results(related_clusters, count_results_by_class, reps_and_alleles_ids,
                                frequency_in_genomes, recommendations, reverse_matches, results_output):
    """
    Writes summary results of BLAST analysis to TSV files.

    This function generates two files: 'related_matches.tsv' and 'count_results_by_cluster.tsv'.
    The 'related_matches.tsv' file contains information about related clusters, while
    'count_results_by_cluster.tsv' details the count of results by cluster and class, including a total count.

    Parameters
    ----------
    related_clusters : dict
        A dictionary where each key is a cluster identifier and each value is a list of tuples.
        Each tuple represents a related match with its details.
    count_results_by_class : dict
        A dictionary where each key is a clusters identifiers separate by '|' and each value is another dictionary.
        The inner dictionary's keys are class identifiers, and values are counts of results for that class.
    reps_and_alleles_ids : dict
        A dictionary mapping pairs of query and subject sequences to their unique loci/CDS IDs and alleles IDs.
    frequency_in_genomes : dict
        A dictionary mapping sequence identifiers to their frequency in genomes.
    reverse_matches : bool
        A flag indicating whether there are reverse matches
    results_output : str
        The path to the directory where the output files will be saved.
    
    Returns
    -------
    None
        This function does not return any value but writes the summary results to the specified files.

    Notes
    -----
    - The 'related_matches.tsv' file is formatted such that each related match is written on a new line,
      with details separated by tabs. A blank line is added after each cluster's matches.
    - The 'count_results_by_cluster.tsv' file includes the cluster identifier, class identifier, count of results,
      and total count of results for the cluster, with each piece of information separated by tabs.
      A blank line is added after each cluster's information.
    """
    related_matches = os.path.join(results_output, "related_matches.tsv")
    reported_cases = {}
    for key, related in list(related_clusters.items()):
        for index, r in enumerate(list(related)):
            if reverse_matches:
                r.insert(4, '-')
                r.insert(5, '-')
            [query, subject] = [itf.remove_by_regex(i, r"\*") for i in r[:2]]
            if (query, subject) not in itf.flatten_list(reported_cases.values()):
                reported_cases.setdefault(key, []).append((subject, query))
            elif reverse_matches:
                sublist_index = itf.find_sublist_index([[itf.remove_by_regex(i, r"\*") for i in l[:2]] for l in related_clusters[key]], [subject, query])
                insert = r[2] if not None else '-'
                related[sublist_index][4] = insert
                insert = r[3] if not None else '-'
                related[sublist_index][5] = insert
                related.remove(r)
        
        for index, i in enumerate(recommendations[key]):
            if index <= (len(related_clusters[key]) - 1):
                related_clusters[key][index] += ([itf.flatten_list([[k] + [i for i in v]]) for k , v in recommendations[key].items()][index])
            else:
                related.append(itf.repeat_strings_in_a_list('', 8)
                               if reverse_matches
                               else itf.repeat_strings_in_a_list('', 6)
                               +
                               ([itf.flatten_list([[k] + [i for i in v]])
                                 for k , v in recommendations[key].items()][index]))

    with open(related_matches, 'w') as related_matches_file:
        related_matches_file.write("Query\tSubject\tClass\tClass_count" +
                                    ("\tInverse_class\tInverse_class_count" if reverse_matches else "") +
                                    "\tFrequency_in_genomes_query\tFrequency_in_genomes_subject\n")
        for related in related_clusters.values():
            for r in related:
                related_matches_file.write('\t'.join(str(item) for item in r) + '\n')

            related_matches_file.write('\n')

    count_results_by_cluster = os.path.join(results_output, "count_results_by_cluster.tsv")
    with open(count_results_by_cluster, 'w') as count_results_by_cluster_file:
        count_results_by_cluster_file.write("Query\tSubject\tClass\tClass_count\tRepresentatives_count"
                                            "\tAlelles_count\tFrequency_in_genomes_query"
                                            "\tFrequency_in_genomes_subject\n")
        for id_, classes in count_results_by_class.items():
            count_results_by_cluster_file.write('\t'.join(id_.split('|')))
            total_count = sum(classes.values())
            query = itf.try_convert_to_type(id_.split('|')[0], int)
            subject = itf.try_convert_to_type(id_.split('|')[1], int)
            for i, items in enumerate(classes.items()):
                if i == 0:
                    count_results_by_cluster_file.write(f"\t{items[0]}\t{items[1]}\{total_count}"
                                                        f"\t{len(reps_and_alleles_ids[id_][0])}"
                                                        f"\t{len(reps_and_alleles_ids[id_][1])}"
                                                        f"\t{frequency_in_genomes[query]}"
                                                        f"\t{frequency_in_genomes[subject]}\n")
                else:
                    count_results_by_cluster_file.write(f"\t\t{items[0]}\t{items[1]}\{total_count}\n")
            count_results_by_cluster_file.write('\n')

def get_loci_matches(all_relationships, loci_ids, cds_to_keep, schema_loci_short, cds_joined_cluster,
                     sorted_blast_dict):
    """
    Determines the matches between loci and their corresponding alleles or CDS based on the
    relationships and the current selection of CDS to keep.

    This function evaluates the relationships between loci and alleles or CDS to identify matches.
    It operates in two modes based on the presence of loci IDs: one where no specific loci IDs are
    provided and it uses all available relationships, and another where specific loci IDs are used
    to filter the matches. It also considers whether the loci or CDS are part of a joined cluster
    and adjusts the matching process accordingly.

    Parameters
    ----------
    all_relationships : dict
        A dictionary containing all relationships between loci and alleles or CDS, with loci as keys
        and lists of related alleles or CDS as values.
    loci_ids : list
        A list of loci IDs to filter the matches. If empty, all relationships are considered.
    cds_to_keep : dict
        A dictionary with classes as keys and lists of CDS or loci IDs to be kept as values.
    schema_loci_short : list
        A list of loci IDs that are considered short and may not be included in the standard matching
        process.
    cds_joined_cluster : dict
        A dictionary mapping integer IDs to their corresponding cluster IDs for CDS that are part of
        a joined cluster.
    sorted_blast_dict : dict
        A dictionary containing sorted BLAST results, used to identify loci that have matches.

    Returns
    -------
    is_matched : dict
        A dictionary with loci or CDS IDs as keys and sets of matched loci IDs as values, indicating
        successful matches.
    is_matched_alleles : dict or None
        A dictionary similar to `is_matched` but specifically for alleles, or None if no specific loci
        IDs are provided.

    Notes
    -----
    - The function first checks if `loci_ids` is provided to determine the mode of operation.
    - It uses utility functions like `itf.flatten_list` to simplify the structure of `all_relationships`
    and `itf.remove_by_regex` to clean up the IDs for matching.
    - The matching process accounts for whether entries are part of a joined cluster and adjusts the
    matching logic accordingly.
    - The function returns two dictionaries: one for general matches and one specifically for alleles,
    the latter being applicable only when `loci_ids` are provided.
    """
    is_matched = {}
    is_matched_alleles = None if not all(loci_ids) else {}
    if not all(loci_ids):
        relationships = itf.flatten_list(all_relationships.values())
        for class_, entries in list(cds_to_keep.items()):
            for entry in list(entries):
                if entry not in schema_loci_short:
                    id_ = entry
                    if isinstance(entry, int):
                        entry = cds_joined_cluster[entry]
                    else:
                        entry = [entry]
                    is_matched.setdefault(id_, set([itf.remove_by_regex(i[0], '_(\d+)') for i in relationships if i[1] in entry]))
    else:
        relationships = itf.flatten_list(all_relationships.values())
        changed_ids = [[r[0], itf.remove_by_regex(r[1], '_(\d+)')] for r in relationships]
        had_matches = set([itf.remove_by_regex(rep, '_(\d+)') for rep in sorted_blast_dict])
        is_matched_alleles = {}
        for class_, entries in list(cds_to_keep.items()):
            for entry in list(entries):
                if entry not in had_matches and not isinstance(entry, int):
                    id_ = entry
                    entry = [entry]
                    is_matched.setdefault(id_, set([i[0] for i in changed_ids if i[1] in entry]))
                    is_matched_alleles.setdefault(id_, set([i[1] 
                                                            for i in relationships 
                                                            if i[0] in is_matched[id_] 
                                                            and itf.remove_by_regex(i[1], '_(\d+)') in entry]))
    return is_matched, is_matched_alleles

def wrap_up_blast_results(cds_to_keep, not_included_cds, clusters, output_path, 
                          constants, drop_set, loci, groups_paths_old, frequency_in_genomes,
                          only_loci):
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
    drop_set : set
        Contains the CDS IDs to be removed from further processing for appearing
        fewer time in genomes than their match.
    loci : dict
        Dict that contains the loci IDs and paths.
    groups_paths_old : dict
        The dictionary containing the old paths for the CDSs groups used 
        to cp instead of creating new FASTAs files.
    frequency_in_genomes : dict
        Dict that contains sum of frequency of that representatives cluster in the
        genomes of the schema.
    only_loci : bool
        If only loci are being processed.

    Returns
    -------
    groups_paths_reps : dict
        Dict that contains as Key the ID of each group while the value is the
        path to the FASTA file that contains its nucleotide sequences.
    reps_trans_dict_cds : dict
        Dict that contais the translations of all the CDSs inside the various
        groups.
    master_file_rep : str or None
        Path to the master file that contains all of the representative sequences.
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
                if only_loci:
                    print(f"\t\tOut of those {count} loci are recommended to be removed.")
                else:
                    print(f"\t\tOut of those {count} {'CDSs groups' if i== 0 else 'loci'}"
                        f" {'were removed from the analysis' if i == 0 else 'are recommended to be replaced with their matched CDS in the schema.'}")
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

        id_folder = os.path.join(output_path, 'results_IDs')
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
            The dictionary containing the old paths for the CDSs groups used 
            to cp instead of creating new FASTAs files.
        loci : dict
            The dictionary containing the loci IDs and paths.

        Returns
        -------
        None, copies file from origin to destination.
        """
        for cds in cds_list:
            if class_ == '1a':
                class_name_cds = f"joined_{cds}"
                for i in cds_list[cds]:
                    file_path = os.path.join(cds_outcome_results, class_name_cds)
                    origin_path = groups_paths_old.pop(i) if case_id == 0 else loci[i]
                    if not os.path.exists(file_path):
                        ff.copy_file(origin_path, file_path)
                    else:
                        ff.concat_files(origin_path, file_path)
                continue
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
            if isinstance(cds,str):
                cds = [cds]
            else:
                cds = cds_to_keep[class_][cds]
            # Write all of the alleles to the file.
            with open(cds_group_fasta_file, 'w') as fasta_file:
                for rep_id in cds:
                    cds_ids = [cds_id for cds_id in clusters[rep_id]]
                    for cds_id in cds_ids:
                        fasta_file.write(f">{cds_id}\n{str(not_included_cds[cds_id])}\n")
            # Write only the representative to the file.
            with open(cds_group_reps_file, 'w') as fasta_file:
                for rep_id in cds:
                    fasta_file.write(f">{rep_id}\n{str(not_included_cds[rep_id])}\n")
            # Write the representative to the master file.
            write_type = 'a' if os.path.exists(master_file_rep) else 'w'
            with open(master_file_rep, write_type) as fasta_file:
                for rep_id in cds:
                    fasta_file.write(f">{rep_id}\n{str(not_included_cds[rep_id])}\n")

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
        reps_trans_dict_cds : dict
            The dictionary containing the translated sequences.
        """
        groups_trans_folder = os.path.join(fasta_folder, 'cds_groups_translation')
        ff.create_directory(groups_trans_folder)
        groups_trans = {}
        for key, group_path in groups_paths.items():
            trans_path = os.path.join(groups_trans_folder, os.path.basename(group_path))
            groups_trans[key] = trans_path
            fasta_dict = sf.fetch_fasta_dict(group_path, False)
            trans_dict, _, _ = sf.translate_seq_deduplicate(fasta_dict,
                                                            trans_path,
                                                            None,
                                                            constants[5],
                                                            False,
                                                            constants[6],
                                                            False)
        
        group_trans_rep_folder = os.path.join(fasta_folder, 'cds_groups_translation_reps')
        ff.create_directory(group_trans_rep_folder)
        groups_trans_reps_paths = {}
        reps_trans_dict_cds = {}
        for key, group_path in groups_paths_reps.items():
            trans_path = os.path.join(group_trans_rep_folder, os.path.basename(group_path))
            groups_trans_reps_paths[key] = trans_path
            fasta_dict = sf.fetch_fasta_dict(group_path, False)
            trans_dict, _, _ = sf.translate_seq_deduplicate(fasta_dict,
                                                            trans_path,
                                                            None,
                                                            constants[5],
                                                            False,
                                                            constants[6],
                                                            False)
            for id_, sequence in trans_dict.items():
                reps_trans_dict_cds[id_] = sequence

        return reps_trans_dict_cds

    def write_cluster_members_to_file(output_path, cds_to_keep, clusters, frequency_in_genomes):
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
        frequency_in_genomes : dict
            Dict that contains sum of frequency of that representatives cluster in the
            genomes of the schema.

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
                        cluster_members_file.write('\t' + str(rep_id))
                        cds_ids = [cds_id for cds_id in clusters[rep_id]]
                        for count, cds_id in enumerate(cds_ids):
                            if count == 0:
                                cluster_members_file.write('\t' + cds_id + '\t' + str(frequency_in_genomes[rep_id]) + '\n')
                            else:
                                cluster_members_file.write('\t\t' + cds_id + '\n')
    # Create directories.
    
    fasta_folder = os.path.join(output_path, 'results_fastas')
    ff.create_directory(fasta_folder)
    
    cds_outcome_results_fastas_folder = os.path.join(fasta_folder, 'results_group_dna_fastas')
    ff.create_directory(cds_outcome_results_fastas_folder)
    if not loci:
        cds_outcome_results_reps_fastas_folder = os.path.join(fasta_folder, 'results_group_dna_reps_fastas')
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

        # Process drop_set in the same way as above
        loci_cases['dropped'] = [d for d in drop_set if d in loci]
        cds_cases['dropped'] = [d for d in drop_set if d not in loci]

    else:
        for class_, cds_set in cds_to_keep.items():
            if class_ == '1a':
                count_cases[class_] = len(itf.flatten_list(cds_set.values()))
            else:
                count_cases[class_] = len(cds_set)

    # Check if loci is not empty
    if loci:
        for i, printout in enumerate([cds_cases, loci_cases]):
            if only_loci and i == 0:
                continue
            total_loci = len(itf.flatten_list([i 
                                               for class_, i
                                               in printout.items()
                                               if class_ != '1a'])) + len(itf.flatten_list(cds_to_keep['1a'].values()))

            print(f"Out of {len(groups_paths_old) if i==0 else len(loci)} {'CDSs groups' if i == 0 else 'loci'}:")
            print(f"\t{total_loci} {'CDSs' if i == 0 else 'loci'}"
                f" representatives had matches with BLASTn against the {'CDSs' if i == 0 else 'schema'}.")

            # Print the classification results
            for class_, group in printout.items():
                print_classification_results(class_ ,len(group) if class_ != '1a' else len(itf.flatten_list(group.values())) ,printout, i)

            if i == 0:
                print(f"\t{len(groups_paths_old) - len(itf.flatten_list(printout.values()))}"
                    " didn't have any BLASTn matches so they were retained.\n")
    else:
        # Write info about the classification results.
        print(f"Out of {len(clusters)} clusters:")
        print(f"\t{sum(count_cases.values()) + len(drop_set)} CDS representatives had matches with BLASTn"
            f" which resulted in {len(itf.flatten_list(cds_to_keep.values()))} groups")

        # Print the classification results
        for class_, count in count_cases.items():
            print_classification_results(class_, count, cds_to_keep, 0)

        print(f"\t\tOut of those {len(drop_set)} CDSs groups were removed from the analysis.")

        if Retained_not_matched_by_blastn:
            print(f"\t{len(Retained_not_matched_by_blastn)} didn't have any BLASTn matches so they were retained.")
            
            cds_to_keep['Retained_not_matched_by_blastn'] = Retained_not_matched_by_blastn
    # Skip the next step to copy or write FASTAS because we are working with the
    # schema only.

    # Initialize dictionaries to store paths
    groups_paths = {}
    groups_paths_reps = {}
    # Check if loci is not None.
    if loci:
        print("\nWriting FASTA file for possible new loci...")
        for case_id, cases in enumerate([cds_cases, loci_cases]):
            if only_loci and case_id == 0:
                continue
            # Create directories and write dict to TSV
            cds_outcome_results = create_directory_and_write_dict(cds_outcome_results_fastas_folder, output_path, case_id, cases)

            # Process each class and CDS list in cases
            for class_, cds_list in cases.items():
                copy_fasta(class_, cds_list, case_id, cds_outcome_results, groups_paths_old, loci)
            # Copy CDS that didnt match
            if case_id == 0 and only_loci:
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
        write_cluster_members_to_file(output_path, cds_to_keep, clusters, frequency_in_genomes)

        master_file_rep = os.path.join(fasta_folder, 'master_rep_file.fasta')

    return groups_paths, reps_trans_dict_cds, master_file_rep

def run_blasts(blast_db, cds_to_blast, reps_translation_dict,
               rep_paths_nuc, output_dir, constants, cpu, multi_fasta = None, if_loci = None):
    """
    This functions runs both BLASTn and Subsequently BLASTp based on results of
    BLASTn.
    
    Parameters
    ----------
    blast_db : str
        Path to the BLAST db folder.
    cds_to_blast : list
        A list of CDS IDs to be used for BLASTn against the BLAST db.
    reps_translation_dict : dict
        A dictionary mapping sequence IDs to their translations (amino acid sequences).
    rep_paths_nuc : dict
        A dictionary mapping CDS IDs to the path of their corresponding FASTA files.
    output_dir : str
        The directory path where output files will be saved.
    constants : list
        A list of constants used within the function, such as thresholds for filtering BLAST results.
    cpu : int
        The number of CPU cores to use for parallel processing.
    multi_fasta : dict, optional
       A dictionary used when the input FASTA files contain multiple CDSs, to ensure correct BLASTn
       execution.
    if_loci : bool, optional
        A flag indicating whether the function should process loci instead of CDSs.
        
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
          "\nRunning BLASTn between groups representatives against schema loci short..." if not if_loci else
          "\nRunning BLASTn between loci representatives against schema loci...")
    # BLASTn folder
    blastn_output = os.path.join(output_dir, '1_BLASTn_processing')
    ff.create_directory(blastn_output)
    # Create directory
    blastn_results_folder = os.path.join(blastn_output, 'BLASTn_results')
    ff.create_directory(blastn_results_folder)
    # Run BLASTn
    # Calculate max id length for print.
    max_id_length = len(max(cds_to_blast))
    total_reps = len(rep_paths_nuc)
    representative_blast_results = {}
    representative_blast_results_coords_all = {}
    representative_blast_results_coords_pident = {}
    # Get Path to the blastn executable
    get_blastn_exec = lf.get_tool_path('blastn')
    i = 1
    with concurrent.futures.ProcessPoolExecutor(max_workers=cpu) as executor:
        for res in executor.map(bf.run_blastdb_multiprocessing,
                                repeat(get_blastn_exec),
                                repeat(blast_db),
                                rep_paths_nuc.values(),
                                cds_to_blast,
                                repeat(blastn_results_folder)
                                ):

            filtered_alignments_dict, _, alignment_coords_all, alignment_coords_pident = af.get_alignments_dict_from_blast_results(
                res[1], constants[1], True, False, True, if_loci)
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
    blastp_results = os.path.join(output_dir, '2_BLASTp_processing')
    ff.create_directory(blastp_results)
    
    blastn_results_matches_translations = os.path.join(blastp_results,
                                                       'blastn_results_matches_translations')
    ff.create_directory(blastn_results_matches_translations)

    representatives_blastp_folder = os.path.join(blastn_results_matches_translations,
                                                'cluster_rep_translation')
    ff.create_directory(representatives_blastp_folder)
    
    blastp_results_folder = os.path.join(blastp_results,
                                         'BLASTp_results')
    ff.create_directory(blastp_results_folder)
    
    blastp_results_ss_folder = os.path.join(blastp_results,
                                            'BLASTp_results_self_score_results')
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
                if isinstance(filename, int):
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
    # Get Path to the blastp executable
    get_blastp_exec = lf.get_tool_path('blastp')
    i = 1
    # Calculate self-score
    with concurrent.futures.ProcessPoolExecutor(max_workers=cpu) as executor:
        for res in executor.map(bf.run_self_score_multiprocessing,
                                rep_paths_prot.keys(),
                                repeat(get_blastp_exec),
                                rep_paths_prot.values(),
                                repeat(blastp_results_ss_folder)):
            
            _, self_score, _, _ = af.get_alignments_dict_from_blast_results(res[1], 0, True, True, True, if_loci)
    
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
        for res in executor.map(bf.run_blast_fastas_multiprocessing,
                                blastp_runs_to_do, 
                                repeat(get_blastp_exec),
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

def write_processed_results_to_file(cds_to_keep, representative_blast_results,
                                    classes_outcome, all_alleles, is_matched,
                                    is_matched_alleles, all_loci, output_path):
    """
    Writes the results of the classification and matching process to files for further analysis or review. 
    This function takes the processed data, including the CDS to keep, representative BLAST results, 
    classification outcomes, allele information, and matching results, and writes them to specified files 
    within a given output directory. It is designed to organize and present the data in a way that facilitates 
    easy access and interpretation.
    
    Parameters
    ----------
    cds_to_keep : dict
        A dictionary categorizing CDS by their classification for retention.
    representative_blast_results : dict
        A nested dictionary with query identifiers as keys, each mapping to another dictionary of subject 
        identifiers and their match details, organized by classification.
    classes_outcome : list
        A list of class IDs, structured as a list of lists, indicating the classification outcomes used 
        in subsequent analyses.
    all_alleles : dict
        A dictionary mapping loci or joined group identifiers to their corresponding alleles and element IDs. 
        This can be `None` if the process does not involve loci.
    is_matched : dict
        A dictionary indicating which CDS/loci have been matched, organized by their identifiers.
    is_matched_alleles : dict
        A dictionary detailing the alleles associated with matched CDS/loci.
    all_loci : bool
        A flag that indicates if only loci are present.
    output_path : str
        The file path to the directory where the output files will be written.

    Creates
    -------
    Files in output directory : 
        Multiple files are created in the specified `output_path`, each containing parts of the processed 
        data. The files are organized to reflect the structure of the input data and the results of the 
        classification and matching process.

    Notes
    -----
    - The function is designed to handle complex data structures resulting from bioinformatics analyses, 
    such as BLAST searches, and to organize this information into a more accessible format.
    - It ensures that the results are not only stored for record-keeping but also formatted in a way that 
    supports easy review and further analysis.
    - The specific format and naming of the output files are determined within the function, based on the 
    structure of the input data and the requirements of the subsequent analysis steps.
    """
    def process_clusters(cds_to_keep, representative_blast_results, all_alleles, is_matched,
                         is_matched_alleles, add_group_column, output_path):
        """
        Processes the results of cluster analysis, specifically focusing on the classification and
        matching of Coding Sequences (CDS) or loci. It iterates through each class of CDS, excluding
        those not matched by BLASTn, to process and document the details of each cluster.
        The function generates a report for each cluster, detailing the CDS or loci involved,
        their match status, and other relevant information. The reports are saved as TSV files
        in a specified output directory. Additionally, the function determines whether an extra
        column for group names is necessary in the report, based on the processed clusters.

        Parameters
        ----------
        cds_to_keep : dict
            A dictionary categorizing CDS by their classification for retention.
        representative_blast_results : dict
            A nested dictionary with query identifiers as keys, each mapping to another dictionary of
            subject identifiers and their match details.
        all_alleles : dict
            A dictionary mapping loci or joined group identifiers to their corresponding alleles and
            element IDs.
        is_matched : dict
            A dictionary indicating which CDS/loci have been matched, organized by their identifiers.
        is_matched_alleles : dict
            A dictionary detailing the alleles associated with matched CDS/loci.
        add_group_column : bool
            A flag indicating whether an additional column for group names should be included in the report.
        output_path : str
            The file path to the directory where the output files will be written.

        Returns
        -------
        None
            Creates and writes TSV files to the specified output directory.

        Notes
        -----
        - The function skips processing for the class 'Retained_not_matched_by_blastn'.
        - It utilizes helper functions such as `process_cluster` to obtain cluster details and
        `generate_write_dict` to prepare data for writing.
        - The output TSV files are named according to the cluster type and its identifier, facilitating
        easy identification and review.
        """
        # Loop over each class and its corresponding CDS
        for class_, cds in cds_to_keep.items():
            if class_ == 'Retained_not_matched_by_blastn':
                continue
            # Loop over each cluster in the CDS
            for id_, cluster in enumerate(cds, 1):
                # Process the cluster and get the necessary details
                id_, cluster, cluster_type, is_cds = process_cluster(class_, id_,
                                                                    cluster,
                                                                    all_alleles,
                                                                    cds)
                # Generate a dictionary to be written to the file
                write_dict = generate_write_dict(id_, cluster, is_cds, is_matched, is_matched_alleles,
                                                 representative_blast_results)
                # Define the path of the report file
                report_file_path = os.path.join(output_path, f"blast_{cluster_type}_{id_}.tsv")
                # Write the dictionary to the file
                alignment_dict_to_file(write_dict, report_file_path, 'w', add_group_column)

    def process_cluster(class_, id_, cluster, all_alleles, cds):
        """
        Processes a single cluster, determining its type, elements, and whether it represents a CDS or a loci. 
        It also identifies if additional group IDs are present, affecting the structure of the output report.

        Parameters
        ----------
        class_ : str
            The classification of the cluster.
        id_ : str or int
            The identifier of the cluster.
        cluster : str or int
            The identifier of the cluster, used when `class_` does not indicate a joined cluster.
        all_alleles : dict
            A dictionary mapping loci or joined group identifiers to their corresponding alleles and
            element IDs.
        cds : dict or str
            Information about the CDS; if it's a single CDS, it contains a string, if it's a joined cluster,
            it contains a dictionary.

        Returns
        -------
        id_ : str or int
            The identifier of the cluster.
        cluster : list
            A list containing the elements' IDs of the cluster.
        cluster_type : str
            The type of the cluster, indicating if it's a joined cluster, retained, CDS cluster, or loci.
        is_cds : bool
            True if the cluster represents a CDS, False if it represents loci.

        Notes
        -----
        - The function first checks the class of the cluster to determine its type and elements.
        - It then assesses whether the cluster represents a CDS or loci based on the presence of alleles
        in `all_alleles`.
        - The presence of additional group IDs is determined by the structure of `all_alleles` and the
        type of entries in the cluster.
        - This function is designed to process clusters in a context where distinguishing between CDS
        and loci, as well as identifying joined clusters, is crucial.
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
            is_cds = False
            cluster_alleles = []
            for entry in cluster:
                if entry not in all_alleles or isinstance(entry, int):
                    if isinstance(entry, int):
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

        return id_, cluster, cluster_type, is_cds

    def generate_write_dict(id_, cluster, is_cds, is_matched, is_matched_alleles,
                            representative_blast_results):
        """
        Generates a dictionary structured for writing to a file, based on the provided cluster
        information, match status, and BLAST results. This function is tailored to handle
        different scenarios, including whether the cluster represents a CDS, if it has been
        matched, and the specifics of those matches.

        Parameters
        ----------
        id_ : str or int
            The identifier of the cluster, which can be a string or an integer.
        cluster : list
            A list containing the identifiers of elements within the cluster.
        is_cds : bool
            A boolean indicating if the cluster represents a Coding Sequence (CDS).
        is_matched : dict
            A dictionary indicating which clusters have been matched, keyed by cluster ID.
        is_matched_alleles : dict
            A dictionary containing the alleles of the matched clusters, keyed by cluster ID.
        representative_blast_results : dict
            A dictionary containing the BLAST results, structured with query identifiers as keys
            and subject identifiers with their match details as values.

        Returns
        -------
        write_dict : dict
            A dictionary formatted for writing to a file. The structure of this dictionary varies
            depending on the match status and type of the cluster (CDS or not).

        Notes
        -----
        - The function handles three main scenarios:
            1. When the cluster represents a CDS and has matches, it compiles a dictionary of these matches,
            filtering by the cluster's elements.
            2. When the cluster itself didn't match but was matched against, it generates a dictionary
            based on the matches and the alleles of the matched clusters.
            3. For all other cases, it creates a dictionary including all subjects for each query within
            the cluster.
        - The function dynamically adjusts the structure of the `write_dict` based on the input parameters,
        ensuring the output is tailored for the specific scenario.
        """
        # Check if it's a CDS and if it matches loci
        if is_cds and is_matched:
            queries = []
            if isinstance(id_, int):
                queries = is_matched[id_]
            else:
                for c in cluster:
                    queries += is_matched[c]
            # Generate the dictionary to be written
            write_dict = {query : {subject: {id_: entry for id_, entry in entries.items()}
                                for subject, entries in subjects.items() if subject in cluster}
                        for query, subjects in representative_blast_results.items()
                        if itf.remove_by_regex(query, '_(\d+)') in queries}
        # for cases that didn't match anything but got matched against.
        elif is_matched and id_ in is_matched:
            queries = is_matched[id_]
            cluster = is_matched_alleles[id_]
            # Generate the dictionary to be written
            write_dict = {query : {subject: {id_: entry for id_, entry in entries.items()}
                                for subject, entries in subjects.items() if subject in cluster}
                        for query, subjects in representative_blast_results.items()
                        if itf.remove_by_regex(query, '_(\d+)') in queries}
        # For all other normal cases.
        else:
            # Generate the dictionary to be written
            write_dict = {query : {subject: {id_: entry for id_, entry in entries.items()}
                                for subject, entries in subjects.items()}
                        for query, subjects in representative_blast_results.items()
                        if query in cluster}
        return write_dict

    def process_classes(classes_outcome, representative_blast_results, output_path, add_group_column):
        """
        Processes the outcomes of different classes from BLAST results and writes the results to files.
        For each class outcome, it generates a dictionary of representative BLAST results filtered by
        class. This dictionary is then written to a TSV file in the specified output directory. The
        function can optionally add a column header for group information based on the `add_group_column`
        parameter.

        Parameters
        ----------
        classes_outcome : list
            A list of class outcomes to process.
        representative_blast_results : dict
            A dictionary containing BLAST results, structured with query identifiers as keys and subject
            identifiers with their match details as values.
        output_path : str
            The path to the directory where the output files will be written.
        add_group_column : bool
            A boolean indicating whether to add a column header for group information in the output files.

        Returns
        -------
        None
            The function does not return any value. It writes the results to TSV files in the specified
            output directory.

        Notes
        -----
        - The function iterates over each class outcome, creating a filtered dictionary of BLAST results
        for that class.
        - It constructs the file path for each class's report using the `output_path` and the class name,
        then writes the filtered results to this file.
        - The `alignment_dict_to_file` function is used to write the dictionary to a TSV file, with the
        option to add a group column if `add_group_column` is True.
        - This function is useful for organizing BLAST results by class and facilitating further analysis
        of these results.
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
            alignment_dict_to_file(write_dict, report_file_path, 'w', add_group_column)

    # Create directories for output
    blast_by_cluster_output = os.path.join(output_path, 'blast_by_cluster')
    ff.create_directory(blast_by_cluster_output)
    blast_results_by_class_output = os.path.join(output_path, 'blast_results_by_class')
    ff.create_directory(blast_results_by_class_output)

    add_group_column = True if not all_loci and all_alleles else False
    # Process and write cluster results
    process_clusters(cds_to_keep, representative_blast_results, all_alleles,
                    is_matched, is_matched_alleles, add_group_column, blast_by_cluster_output)

    # Process and write class results
    process_classes(classes_outcome, representative_blast_results, blast_results_by_class_output,
                    add_group_column)

def extract_cds_to_keep(classes_outcome, count_results_by_class, drop_mark):
    """
    Extracts and organizes CDS (Coding Sequences) to keep based on classification outcomes.

    This function processes BLAST results to determine which coding sequences (CDS) should
    be retained for further analysis based on their classification outcomes. It organizes
    CDS into categories, prioritizes them according to a predefined order of classes, and
    identifies sequences to be dropped.

    Parameters
    ----------
    classes_outcome : list
        An ordered list of class identifiers that determine the priority of classes for keeping CDS.
    count_results_by_class : dict
        A dictionary where keys are concatenated query and subject IDs separated by '|', and values
        are dictionaries with class identifiers as keys and counts as values.
    drop_mark : set
        A set of identifiers that are marked for dropping based on previous criteria.
    cds_to_keep : dict
        An initially empty dictionary that will be populated with CDS to keep, organized by class.

    Returns
    -------
    cds_to_keep : dict
        A dictionary with class identifiers as keys and lists of CDS identifiers or pairs of identifiers
        to be kept in each class.
    drop_set : set
        A set of CDS identifiers that are determined to be dropped based on their classification and
        presence in `drop_mark`.

    Notes
    -----
    - The function first initializes `cds_to_keep` with empty lists for each class in `classes_outcome`.
    - It then iterates through `count_results_by_class` to assign CDS to the most appropriate class
    based on the provided outcomes.
    - Special handling is given to class '1a', where CDS pairs are clustered and indexed.
    - CDS marked in `drop_mark` and falling under certain classes are added to `drop_set` for exclusion.
    - The function uses utility functions like `itf.try_convert_to_type` for type conversion and
    `cf.cluster_by_ids` for clustering CDS pairs in class '1a'.
    """
    temp_keep = {}
    cds_to_keep = {}
    drop_set = set()
    for class_ in classes_outcome:
        cds_to_keep[class_] = []
    for ids, result in count_results_by_class.items():
        class_ = next(iter(result))
        [query, subject] = list(map(lambda x: itf.try_convert_to_type(x, int), ids.split('|')))
        if class_ == '1a':
            cds_to_keep.setdefault('1a', []).append([query, subject])
        if not temp_keep.get(query):
            temp_keep[query] = class_
        elif classes_outcome.index(class_) < classes_outcome.index(temp_keep[query]):
            temp_keep[query] = class_
        if not temp_keep.get(subject):
            temp_keep[subject] = class_
        elif classes_outcome.index(class_) < classes_outcome.index(temp_keep[subject]):
            temp_keep[subject] = class_

    for keep, class_ in temp_keep.items():
        if class_ == '1a':
            continue
        if keep in drop_mark and class_ in ['1b', '2a', '3a']:
            drop_set.add(itf.try_convert_to_type(keep, int))
        else:
            cds_to_keep.setdefault(class_, []).append(keep)

    cds_to_keep['1a'] = {i: values for i, values in enumerate(cf.cluster_by_ids(cds_to_keep['1a']), 1)}

    return cds_to_keep, drop_set

def process_schema(schema, groups_paths, results_output, reps_trans_dict_cds, 
                   cds_to_keep, frequency_in_genomes, allelecall_directory, 
                   master_file_rep, loci_ids, master_alleles, constants, cpu):
    """
    This function processes data related to the schema seed, importing, translating
    and BLASTing against the unclassified CDS clusters representatives groups to
    validate them.
    
    Parameters
    ----------
    schema : str
        Path to the schema seed folder.
    groups_paths : dict
        Dict that contains the path to the FASTA file for each group.
    results_output : str
        Path were to write the results of this function.
    reps_trans_dict_cds : dict
        Dict that contains the translations for each CDS.
    cds_to_keep : dict     
        Dict of the CDS to keep by each classification.
    frequency_in_genomes : dict
        Dict that contains sum of frequency of that representatives cluster in the
        genomes of the schema.
    allelecall_directory : str
        Path to the allele call directory.
    master_file_rep : str
        Path to the master file containing retained CDS.
    loci_ids : list
        List containg two bools, each representing query and subject, True
        if they are loci False otherwise.
    master_alleles : bool
        If True, the function will process all of the alleles of the loci, if False only the
        representatives.
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
    blast_results = os.path.join(results_output, '1_BLAST_processing')
    ff.create_directory(blast_results)
    # Create BLASTn_processing directory
    blastn_output = os.path.join(blast_results, '1_BLASTn_processing')
    ff.create_directory(blastn_output)

    # Get all of the schema loci short FASTA files path.
    schema_short_path = os.path.join(schema, 'short')
    schema_loci_short = {loci_path.replace("_short.fasta", ""): os.path.join(schema_short_path, loci_path) 
                         for loci_path in os.listdir(schema_short_path) 
                         if loci_path.endswith('.fasta')}
    
    # Get all of the schema loci FASTA files path.
    schema_loci = {loci_path.replace(".fasta", ""): os.path.join(schema, loci_path) 
                         for loci_path in os.listdir(schema) 
                         if loci_path.endswith('.fasta')}

    # Create a folder for short translations.
    blastp_output =  os.path.join(blast_results, '2_BLASTp_processing')
    ff.create_directory(blastp_output)
    short_translation_folder = os.path.join(blastp_output, 'short_translation_folder')
    ff.create_directory(short_translation_folder)

    # Find the file in the allele call results that contains the total of each.
    # classification obtained for each loci.
    results_statistics = os.path.join(allelecall_directory, 'loci_summary_stats.tsv')
    # Convert TSV table to dict.
    results_statistics_dict = itf.tsv_to_dict(results_statistics)
    # Add the results for all of the Exact matches to the frequency_in_genomes dict.
    for key, value in results_statistics_dict.items():
        frequency_in_genomes.setdefault(key, int(value[0]))
    # Translate each short loci and write to master fasta.
    i = 1
    len_short_folder = len(schema_loci_short)
    all_alleles = {}
    if not master_file_rep:
        filename = 'master_file' if master_alleles else 'master_file_rep'
        master_file_rep_folder = os.path.join(blastn_output, filename)
        ff.create_directory(master_file_rep_folder)
        master_file_rep = os.path.join(master_file_rep_folder, f"{filename}.fasta")
        write_to_master = True
    else:
        write_to_master = False
    # Create varible to store proteins sequences if it doesn't exist.
    reps_trans_dict_cds = {} if not reps_trans_dict_cds else reps_trans_dict_cds
    # If to BLAST against reps or all of the alleles.
    schema_loci if master_alleles else schema_loci_short
    for loci, loci_short_path in schema_loci.items():
        print(f"\rTranslated{'' if master_alleles else ' short'} loci FASTA: {i}/{len_short_folder}", end='', flush=True)
        i += 1
        fasta_dict = sf.fetch_fasta_dict(loci_short_path, False)
        
        for allele_id, sequence in fasta_dict.items():
            all_alleles.setdefault(loci, []).append(allele_id)

            if write_to_master:
                write_type = 'w' if not os.path.exists(master_file_rep) else 'a'
                with open(master_file_rep, write_type) as master_file:
                    master_file.write(f">{allele_id}\n{sequence}\n")

        loci_short_translation_path = os.path.join(short_translation_folder, f"{loci}.fasta")
        translation_dict, _, _ = sf.translate_seq_deduplicate(fasta_dict, 
                                                              loci_short_translation_path,
                                                              None,
                                                              constants[5],
                                                              False,
                                                              constants[6],
                                                              False)
        for allele_id, sequence in translation_dict.items():
            reps_trans_dict_cds[allele_id] = sequence

    # Create BLAST db for the schema DNA sequences.
    print(f"\nCreate BLAST db for the {'schema' if master_alleles else 'unclassified'} DNA sequences...")
    makeblastdb_exec = lf.get_tool_path('makeblastdb')
    blast_db = os.path.join(blastn_output, 'blast_db_nucl')
    ff.create_directory(blast_db)
    blast_db_nuc = os.path.join(blast_db, 'Blast_db_nucleotide')
    bf.make_blast_db(makeblastdb_exec, master_file_rep, blast_db_nuc, 'nucl')

    [representative_blast_results,
     representative_blast_results_coords_all,
     representative_blast_results_coords_pident,
     bsr_values,
     _] = run_blasts(blast_db_nuc,
                     schema_loci_short,
                     reps_trans_dict_cds,
                     schema_loci_short,
                     blast_results,
                     constants,
                     cpu,
                     all_alleles,
                     all(loci_ids))

    add_items_to_results(representative_blast_results,
                         None,
                         bsr_values,
                         representative_blast_results_coords_all,
                         representative_blast_results_coords_pident,
                         frequency_in_genomes,
                         loci_ids,
                         cds_to_keep['1a'] if cds_to_keep else None)

    # Add CDS joined clusters to all_alleles IDS
    cds_joined_cluster = cds_to_keep['1a'] if cds_to_keep else None
    if cds_joined_cluster:
        all_alleles.update(cds_joined_cluster)
    # Separate results into different classes.
    classes_outcome = separate_blastn_results_into_classes(representative_blast_results,
                                                           constants)
    report_file_path = os.path.join(results_output, 'blast_all_matches.tsv')
    # Write all of the BLASTn results to a file.
    alignment_dict_to_file(representative_blast_results, report_file_path, 'w', True)
    
    print("\nProcessing classes...")
    sorted_blast_dict = sort_blast_results_by_classes(representative_blast_results, classes_outcome)
    # Process the results_outcome dict and write individual classes to TSV file.
    processed_results, count_results_by_class, reps_and_alleles_ids, drop_mark = process_classes(sorted_blast_dict,
                                                                                        classes_outcome,
                                                                                        all_alleles)
    # Sort the count_results_by_class dict by the classes_outcome tuple.
    count_results_by_class = itf.sort_subdict_by_tuple(count_results_by_class, classes_outcome)
    # Extract CDS to keep and drop set.
    cds_to_keep, drop_set = extract_cds_to_keep(classes_outcome, count_results_by_class, drop_mark)
    # Extract the related clusters and recommendations what to do with them.
    all_relationships, related_clusters, recommendations  = extract_results(processed_results,
                                                                           count_results_by_class,
                                                                           frequency_in_genomes,
                                                                           cds_to_keep,
                                                                           drop_set,
                                                                           all(loci_ids),
                                                                           classes_outcome)

    write_blast_summary_results(related_clusters,
                                count_results_by_class,
                                reps_and_alleles_ids,
                                frequency_in_genomes,
                                recommendations,
                                all(loci_ids),
                                results_output)

    # Get all of the CDS that matched with loci
    [is_matched, is_matched_alleles] = get_loci_matches(all_relationships,
                                                        loci_ids,
                                                        cds_to_keep,
                                                        schema_loci_short,
                                                        cds_joined_cluster,
                                                        sorted_blast_dict)

    print("\nWritting classes results to files...")
    write_processed_results_to_file(cds_to_keep,
                                    sorted_blast_dict,
                                    classes_outcome,
                                    all_alleles,
                                    is_matched,
                                    is_matched_alleles,
                                    all(loci_ids),
                                    results_output)
    
    print("\nWrapping up BLAST results...")

    wrap_up_blast_results(cds_to_keep,
                        None,
                        all_alleles,
                        results_output,
                        constants,
                        drop_set,
                        schema_loci,
                        groups_paths,
                        frequency_in_genomes,
                        all(loci_ids))

    return sorted_blast_dict

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

def classify_cds(schema, output_directory, allelecall_directory, constants, temp_paths, cpu):

    temp_folder = temp_paths[0]
    file_path_cds = temp_paths[1]
    #missing_classes_fastas = temp_paths[2]

    # Verify if the dataset is small, if it is, keep minimum genomes in which
    # specific CDS cluster is present to 5 if not to 1% of the dataset size.
    if not constants[2]:
        count_genomes_path = os.path.join(temp_folder, '1_cds_prediction')
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

    """
    print("Identifying CDS identified as missing classes...")
    missing_classes_fastas = sf.fetch_fasta_dict(missing_classes_fastas, True)
    print("Filtering missing CDS in the schema...")
    missing_classes_fastas = {itf.remove_by_regex(key.split('|')[3], '&.*'): value 
                              for key, value in missing_classes_fastas.items() 
                              if itf.regex_present(['&ASM', '&ALM', '&NIPH', '&NIPHEM'], key)}
    # Deduplicate the FASTA dict.
    missing_classes_fastas = itf.deduplicate_fasta_dict(missing_classes_fastas)

    not_included_cds.update(missing_classes_fastas)
    """
    print("Filtering missing CDS in the schema...")
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

    cds_output = os.path.join(output_directory, '1_CDS_processing')
    ff.create_directory(cds_output)
    # This file contains unique CDS.
    cds_not_present_file_path = os.path.join(cds_output, 'CDS_not_found.fasta')
    
    # Count the number of CDS present in the schema and write CDS sequence
    # into a FASTA file.
    frequency_cds = {}
    with open(cds_not_present_file_path, 'w+') as cds_not_found:
        for id_, sequence in not_included_cds.items():
            cds_not_found.write(f">{id_}\n{str(sequence)}\n")
            
            hashed_seq = sf.seq_to_hash(str(sequence))
            # if CDS sequence is present in the schema count the number of
            # genomes that it is found minus 1 (subtract the first CDS genome).
            if hashed_seq in decoded_sequences_ids:
                frequency_cds[id_] = len(decoded_sequences_ids[hashed_seq]) - 1
            else:
                frequency_cds[id_] = 0
                

    print("\nTranslate and deduplicate CDS...")
    # Translate the CDS and find unique proteins using hashes, the CDS with
    # the same hash will be added under that hash in protein_hashes.
    cds_not_present_trans_file_path = os.path.join(cds_output, "CDS_not_found_translation.fasta")
    cds_not_present_untrans_file_path = os.path.join(cds_output, "CDS_not_found_untranslated.fasta")
    # Translate and deduplicate protein sequences.
    cds_translation_dict, protein_hashes, _ = sf.translate_seq_deduplicate(not_included_cds,
                                                                           cds_not_present_trans_file_path,
                                                                           cds_not_present_untrans_file_path,
                                                                           constants[5],
                                                                           True,
                                                                           constants[6],
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
    frequency_in_genomes = {rep: sum([frequency_cds[entry] for entry in value]) 
                             for rep, value in clusters.items()}
    # Filter cluster by the total sum of CDS that are present in the genomes, based on input value.
    clusters = {rep: cluster_member for rep, cluster_member in clusters.items() 
                if frequency_in_genomes[rep] >= constants[2]}
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
    blast_output = os.path.join(output_directory, '2_BLAST_processing')
    ff.create_directory(blast_output)
    
    blastn_output = os.path.join(blast_output, '1_BLASTn_processing')
    ff.create_directory(blastn_output)
    # Create directory and files path where to write FASTAs.
    representatives_blastn_folder = os.path.join(blastn_output,
                                                'cluster_representatives_fastas')
    ff.create_directory(representatives_blastn_folder)

    representatives_all_fasta_file = os.path.join(representatives_blastn_folder,
                                                  'all_cluster_representatives.fasta')
    # Write files for BLASTn.
    rep_paths_nuc = {}
    # Write master file for the representatives.
    with open(representatives_all_fasta_file, 'w') as all_fasta:
        for cluster_rep_id in clusters:
            all_fasta.write(f">{cluster_rep_id}\n{str(not_included_cds[cluster_rep_id])}\n")

            rep_fasta_file = os.path.join(representatives_blastn_folder,
                                          f"cluster_rep_{cluster_rep_id}.fasta")
            rep_paths_nuc[cluster_rep_id] = rep_fasta_file
            # Write the representative FASTA file.
            with open(rep_fasta_file, 'w') as rep_fasta:
                rep_fasta.write(f">{cluster_rep_id}\n{str(not_included_cds[cluster_rep_id])}\n")
    
    # Create BLAST db for the schema DNA sequences.
    print("\nCreating BLASTn database for the unclassified and missed CDSs...")
    # Get the path to the makeblastdb executable.
    makeblastdb_exec = lf.get_tool_path('makeblastdb')
    blast_db = os.path.join(blastn_output, 'blast_db_nucl', 'blast_nucleotide_db')
    bf.make_blast_db(makeblastdb_exec, representatives_all_fasta_file, blast_db, 'nucl')

    # Run the BLASTn and BLASTp
    [representative_blast_results,
     representative_blast_results_coords_all,
     representative_blast_results_coords_pident,
     bsr_values,
     _] = run_blasts(blast_db,
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
                         frequency_in_genomes,
                         [False, False])

    print("\nFiltering BLAST results into classes...")
    results_output = os.path.join(output_directory, '3_CDS_processing_results')
    ff.create_directory(results_output)
    report_file_path = os.path.join(results_output, 'blast_all_matches.tsv')
    
    # Separate results into different classes.
    classes_outcome = separate_blastn_results_into_classes(representative_blast_results,
                                                           constants)
    # Write all of the BLASTn results to a file.
    alignment_dict_to_file(representative_blast_results, report_file_path, 'w')
    
    print("Processing classes...")
    # Process the results_outcome dict and write individual classes to TSV file.
    processed_results, count_results_by_class, reps_and_alleles_ids, drop_mark = process_classes(representative_blast_results,
                                                                                classes_outcome,
                                                                                None)

    count_results_by_class = itf.sort_subdict_by_tuple(count_results_by_class, classes_outcome)

    cds_to_keep, drop_set = extract_cds_to_keep(classes_outcome, count_results_by_class, drop_mark)

    all_relationships, related_clusters, recommendations = extract_results(processed_results,
                                                                          count_results_by_class,
                                                                          frequency_in_genomes,
                                                                          cds_to_keep,
                                                                          drop_set,
                                                                          True,
                                                                          classes_outcome)
    
    write_blast_summary_results(related_clusters,
                                count_results_by_class,
                                reps_and_alleles_ids,
                                frequency_in_genomes,
                                recommendations,
                                True,
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
                                    None,
                                    [False, False],
                                    results_output)
    
    print("\nWrapping up BLAST results...")
    [groups_paths,
     reps_trans_dict_cds,
     master_file_rep] = wrap_up_blast_results(cds_to_keep,
                                              not_included_cds,
                                              clusters,
                                              results_output,
                                              constants,
                                              drop_set,
                                              None,
                                              None,
                                              frequency_in_genomes,
                                              False)
    
    # Add new frequencies in genomes for joined groups
    new_cluster_freq = {}
    for cluster_id, cluster_members in cds_to_keep['1a'].items():
        new_cluster_freq[cluster_id] = 0
        for member in cluster_members:
            new_cluster_freq[(cluster_id)] += frequency_in_genomes[member]
        for member in cluster_members:
            frequency_in_genomes[member] = new_cluster_freq[cluster_id]
    frequency_in_genomes.update(new_cluster_freq)

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
    results_output = os.path.join(output_directory, '4_Schema_processing')
    ff.create_directory(results_output)

    loci_ids = [True, False]
    # Run Blasts for the found loci against schema short
    representative_blast_results = process_schema(schema,
                                                  groups_paths,
                                                  results_output,
                                                  reps_trans_dict_cds,
                                                  cds_to_keep,
                                                  frequency_in_genomes,
                                                  allelecall_directory, 
                                                  master_file_rep,
                                                  loci_ids,
                                                  False,
                                                  constants,
                                                  cpu)