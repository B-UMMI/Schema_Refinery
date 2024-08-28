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
                         frequency_in_genomes, allele_ids, add_groups_ids):
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
    allele_ids : list
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
        -----add_items_to_results
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
                       local_palign_min, allele_ids, add_groups_ids):
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
        allele_ids : list
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
        if allele_ids[0]:
            query_before = query
            query = itf.remove_by_regex(query, '_(\d+)')
        if allele_ids[1]:
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
        if allele_ids[0]:
            query = query_before
        if allele_ids[1]:
            subject = subject_before
        representative_blast_results[query][subject][entry_id].update(update_dict)

        if add_groups_ids:
            id_ = itf.identify_string_in_dict_get_key(subject, add_groups_ids) or subject
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
                    update_results(representative_blast_results, query, subject, entry_id, bsr, sim, cov, frequency_in_genomes, global_palign_all_min, global_palign_all_max, global_palign_pident_min, global_palign_pident_max, local_palign_min, allele_ids, add_groups_ids)
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
    pident = constants[1]
    bsr = constants[7]
    size_ratio =  1 - constants[9]
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
                if blastn_entry['global_palign_all_min'] >= size_ratio:
                    if blastn_entry['bsr'] >= bsr:
                        # Add to class '1a' if bsr is greater than or equal to bsr value
                        add_class_to_dict('1a')
                    elif freq_ratio <= 0.1:
                        # Add to class '1b' if frequency ratio is less than or equal to 0.1
                        add_class_to_dict('1b')
                    else:
                        # Add to class '1c' if none of the above conditions are met
                        add_class_to_dict('1c')
                elif 0.4 <= blastn_entry['global_palign_all_min'] < size_ratio:
                    if blastn_entry['pident'] >= pident:
                        if blastn_entry['global_palign_pident_max'] >= size_ratio:
                            # Add to class '2a' or '2b' based on frequency ratio
                            add_class_to_dict('2a' if freq_ratio <= 0.1 else '2b')
                        else:
                            # Add to class '3a' or '3b' based on frequency ratio
                            add_class_to_dict('3a' if freq_ratio <= 0.1 else '3b')
                    else:
                        if blastn_entry['global_palign_pident_max'] >= size_ratio:
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
    count_results_by_class_with_inverse : dict
        A dictionary containing counts of results by class, including inverse matches, with keys formatted as
        "query|subject" and values being dictionaries with class identifiers as keys and counts as values.
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
    count_results_by_class_with_inverse = {}
    reps_and_alleles_ids = {}
    processed_results = {}
    all_relationships = {class_: [] for class_ in classes_outcome}
    drop_mark = []
    inverse_match = []
    # Process the CDS to find what CDS to retain while also adding the relationships between different CDS
    for query, rep_blast_result in representative_blast_results.items():
        for id_subject, matches in rep_blast_result.items():
            class_ = matches[1]['class']
            all_relationships[class_].append([query, id_subject])
            ids_for_relationship = [query, id_subject]
            new_query = query
            new_id_subject = id_subject

            strings = [str(query), str(id_subject), class_]
            if all_alleles:
                replaced_query = itf.identify_string_in_dict_get_key(query, all_alleles)
                if replaced_query:
                    new_query = replaced_query
                    strings[0] = new_query
                replaced_id_subject = itf.identify_string_in_dict_get_key(id_subject, all_alleles)
                if replaced_id_subject:
                    new_id_subject = replaced_id_subject
                    strings[1] = new_id_subject

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
            
            if f"{new_query}|{new_id_subject}" not in inverse_match:
                count_results_by_class_with_inverse.setdefault(f"{new_query}|{new_id_subject}", {})
                inverse_match.append(f"{new_id_subject}|{new_query}")
            if f"{new_query}|{new_id_subject}" in inverse_match:
                if not count_results_by_class_with_inverse[f"{new_id_subject}|{new_query}"].get(class_):
                    count_results_by_class_with_inverse[f"{new_id_subject}|{new_query}"].setdefault(class_, ['-', 1])
                elif count_results_by_class_with_inverse[f"{new_id_subject}|{new_query}"][class_][1] == '-':
                    count_results_by_class_with_inverse[f"{new_id_subject}|{new_query}"][class_][1] = 1
                else:
                    count_results_by_class_with_inverse[f"{new_id_subject}|{new_query}"][class_][1] += 1
            else:
                if not count_results_by_class_with_inverse[f"{new_query}|{new_id_subject}"].get(class_):
                    count_results_by_class_with_inverse[f"{new_query}|{new_id_subject}"].setdefault(class_, [1, '-'])
                elif count_results_by_class_with_inverse[f"{new_query}|{new_id_subject}"][class_][0] == '-':
                    count_results_by_class_with_inverse[f"{new_query}|{new_id_subject}"][class_][0] = 1
                else:
                    count_results_by_class_with_inverse[f"{new_query}|{new_id_subject}"][class_][0] += 1
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

    return processed_results, count_results_by_class, count_results_by_class_with_inverse, reps_and_alleles_ids, drop_mark, all_relationships

def extract_results(processed_results, count_results_by_class, frequency_in_genomes,
                    clusters_to_keep, drop_possible_loci, classes_outcome):
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
        A list of class outcomes.process_id

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
    - It uses helper functions like `cf.cluster_by_ids` for clustering and `itf.identify_string_in_dict_get_key`
    for identifying if a query or subject ID is present in the clusters.
    """
    def cluster_data(processed_results):
        """
        Cluster data based on a specific key extracted from the processed results.

        Parameters
        ----------
        processed_results : dict
            A dictionary containing the processed results to be clustered.

        Returns
        -------
        return : dict
            A dictionary where each key is an integer starting from 1, and each value is a cluster of data
            based on the extracted key, excluding entries with specific identifiers.
        """
        key_extractor = lambda v: v[3]
        condition = lambda v: v[0] not in ['4c', '5']

        return {i: cluster for i, cluster in enumerate(cf.cluster_by_ids([key_extractor(v) for v in processed_results.values() if condition(v)]), 1)}

    def choice_data(processed_results, to_cluster_list):
        """
        Select and cluster data based on specific conditions and a list of identifiers to cluster.

        Parameters
        ----------
        processed_results : dict
            A dictionary containing the processed results to be selected and clustered.
        to_cluster_list : list
            A list of identifiers indicating which entries should be considered for clustering.

        Returns
        -------
        return : dict
            A dictionary where each key is an integer starting from 1, and each value is a cluster of data
            selected based on the given conditions and the to_cluster_list.
        """
        key_extractor = lambda v: v[3]
        additional_condition = lambda v: '*' in v[4][0] or '*' in v[4][1]
        return {i: cluster for i, cluster in enumerate(cf.cluster_by_ids([key_extractor(v) for v in processed_results.values() if v[0] not in ['1a','4c','5'] or ((itf.identify_string_in_dict_get_key(v[3][0], to_cluster_list) and additional_condition(v)) or (itf.identify_string_in_dict_get_key(v[3][1], to_cluster_list) and additional_condition(v)))]), 1)}
    
    def process_id(id_, to_cluster_list, clusters_to_keep):
        """
        Process an identifier to check its presence in specific lists.

        Parameters
        ----------
        id_ : str
            The identifier to be processed.
        to_cluster_list : list
            A list of identifiers to check for the presence of id_.
        clusters_to_keep : dict
            A dictionary containing identifiers to check for a joined condition.

        Returns
        -------
        return : tuple
            A tuple containing the original id, a boolean indicating if the id is present in the to_cluster_list,
            and a boolean indicating if the id is present in the clusters_to_keep under a specific key.
        """
        present = itf.identify_string_in_dict_get_key(id_, to_cluster_list)
        joined_id = itf.identify_string_in_dict_get_key(id_, clusters_to_keep['1a'])
        return id_, present, joined_id

    def check_in_recommendations(id_, joined_id, recommendations, key, categories):
        """
        Check if an identifier or its joined form is present in a set of recommendations.

        Parameters
        ----------
        id_ : str
            The original identifier to check.
        joined_id : str
            The joined form of the identifier to check.
        recommendations : dict
            A dictionary of recommendations to search within.
        key : str
            The key within the recommendations to search under.
        categories : list
            A list of categories to consider when searching in the recommendations.

        Returns
        -------
        return : bool
            True if the identifier or its joined form is present in the recommendations under the specified categories,
            False otherwise.
        """
        return any((joined_id or id_) in itf.flatten_list([v for k, v in recommendations[key].items() if cat in k]) for cat in categories)


    def add_to_recommendations(category, id_to_write, joined_id=None):
        """
        Add an identifier to the recommendations under a specific category.

        Parameters
        ----------
        category : str
            The category under which to add the identifier.
        id_to_write : str
            The identifier to add to the recommendations.
        joined_id : str, optional
            The joined form of the identifier, if applicable. Default is None.

        Returns
        -------
        None
            This function does not return any value but modifies the `recommendations` dictionary in place.
        """
        if joined_id is not None:  # For joined or choice categories
            recommendations[key].setdefault(f'{category}_{joined_id}', set()).add(id_to_write)
        else:  # For keep or drop categories
            recommendations[key].setdefault(category, set()).add(id_to_write)

    related_clusters = {} # To keep track of the related clusters
    recommendations = {} # To keep track of the recommendations
    dropped_match = {} # To keep track of the dropped matches
    matched_with_dropped = {} # To keep track of the matches that matched with dropped
    processed_cases = [] # To keep track of the processed cases
    # Normal run, where IDs are only loci or CDS original IDs.
    to_cluster_list = cluster_data(processed_results) # All of the cluster of related CDS/loci by ID.
    choice = choice_data(processed_results, to_cluster_list) # All of the possible cases of choice.

    related_clusters = {}
    for results in processed_results.values():
        if results[0] in ['4c','5', 'Retained_not_matched_by_blastn']:
            continue

        query_id, query_present, joined_query_id = process_id(results[3][0], to_cluster_list, clusters_to_keep)
        subject_id, subject_present, joined_subject_id = process_id(results[3][1], to_cluster_list, clusters_to_keep)

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
    
        if_query_dropped = (joined_query_id or query_id) in drop_possible_loci
        if_subject_dropped = (joined_subject_id or subject_id) in drop_possible_loci

        choice_query_id = itf.identify_string_in_dict_get_key(query_id, choice)
        choice_subject_id = itf.identify_string_in_dict_get_key(subject_id, choice)

        # What IDs to addto the Keep, Drop and Choice.
        query_to_write = joined_query_id or query_id
        subject_to_write = joined_subject_id or subject_id
        
        query_to_write = query_to_write
        subject_to_write = subject_to_write

        joined_query_to_write = query_id
        joined_subject_to_write = subject_id

        # Check if the pair was not processed yet
        if [query_id, subject_id] not in processed_cases:
            processed_cases.append([subject_id, query_id])  # Add the inverse pair to the processed cases
            reverse_id = f"{subject_id}|{query_id}"
            reverse_results = processed_results[reverse_id] if processed_results.get(reverse_id) else None
            if reverse_results and classes_outcome.index(results[0]) > classes_outcome.index(reverse_results[0]):
                results = reverse_results

            if results[0] == '1a':
                if joined_query_id is not None:
                    add_to_recommendations('Joined', joined_query_to_write, joined_query_id)
                if joined_subject_id is not None:
                    add_to_recommendations('Joined', joined_subject_to_write, joined_subject_id)

            elif results[0] in ['1c', '2b', '3b', '4b']:
                if not if_query_dropped and not if_subject_dropped and not if_same_joined:
                    add_to_recommendations('Choice', query_to_write, choice_query_id)
                    add_to_recommendations('Choice', subject_to_write, choice_subject_id)
                elif if_query_dropped:
                    add_to_recommendations('Choice', subject_to_write, choice_subject_id)
                    matched_with_dropped.setdefault(key, []).append([query_to_write, subject_to_write, query_to_write, choice_query_id])
                elif if_subject_dropped:
                    add_to_recommendations('Choice', query_to_write, choice_query_id)
                    matched_with_dropped.setdefault(key, []).append([query_to_write, subject_to_write, subject_to_write, choice_query_id])
            elif results[0] in ['1b', '2a', '3a', '4a']:
                if (joined_query_id and '*' in results[4][0]) or (joined_subject_id and '*' in results[4][1]):
                    add_to_recommendations('Choice', query_to_write, choice_query_id)
                    add_to_recommendations('Choice', subject_to_write, choice_subject_id)
                if query_id in drop_possible_loci:
                    if not if_joined_query and not if_query_in_choice:
                        add_to_recommendations('Drop', query_to_write)
                        dropped_match.setdefault(key, []).append([query_to_write, subject_to_write, subject_to_write])

                elif subject_id in drop_possible_loci:
                    if not if_joined_subject and not if_subject_in_choice:
                        add_to_recommendations('Drop', subject_to_write)
                        dropped_match.setdefault(key, []).append([subject_to_write, query_to_write, query_to_write])
    # Add cases where some ID matched with dropped ID, we need to add the ID that matched with the ID that made the other match
    # to be Dropped. e.g x and y matched with x dropping and x also matched with z, then we need to make a choice between x and z.
    for key, matches in matched_with_dropped.items():
        if dropped_match.get(key):
            for dropped in dropped_match[key]:
                for match_ in  matches:
                    if match_[2] in dropped:
                        add_to_recommendations('Choice', dropped[2], match_[3])
            

    sort_order = ['Joined', 'Choice', 'Keep', 'Drop']
    recommendations = {k: {l[0]: l[1] for l in sorted(v.items(), key=lambda x: sort_order.index(x[0].split('_')[0]))} for k, v in recommendations.items()}
    
    return related_clusters, recommendations

def write_blast_summary_results(related_clusters, count_results_by_class, group_reps_ids, group_alleles_ids,
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
    recommendations : dict
        A dictionary containing recommendations for each cluster based on the classification of the results.
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
            related_clusters[key][index] += ([itf.flatten_list([[k] + [i for i in v]]) for k , v in recommendations[key].items()][index])
    # Write the results to the output files
    related_matches = os.path.join(results_output, "related_matches.tsv")
    with open(related_matches, 'w') as related_matches_file:
        related_matches_file.write("Query\tSubject\tClass\tClass_count" +
                                    ("\tInverse_class\tInverse_class_count" if reverse_matches else "") +
                                    "\tFrequency_in_genomes_query\tFrequency_in_genomes_subject\n")
        for related in related_clusters.values():
            for r in related:
                related_matches_file.write('\t'.join(str(item) for item in r) + '\n')

            related_matches_file.write('#\n')

    count_results_by_cluster = os.path.join(results_output, "count_results_by_cluster.tsv")
    with open(count_results_by_cluster, 'w') as count_results_by_cluster_file:
        count_results_by_cluster_file.write("Query\tSubject\tClass\tClass_count\tInverse_class_count"
                                            "\tRepresentatives_count"
                                            "\tAlelles_count\tFrequency_in_genomes_query"
                                            "\tFrequency_in_genomes_subject\n")
        for id_, classes in count_results_by_class.items():
            query, subject = id_.split('|')
            count_results_by_cluster_file.write('\t'.join(id_.split('|')))
            total_count_origin = sum([i[0] for i in classes.values() if i[0] != '-'])
            total_count_inverse = sum([i[1] for i in classes.values() if i[1] != '-'])
            query = itf.try_convert_to_type(id_.split('|')[0], int)
            subject = itf.try_convert_to_type(id_.split('|')[1], int)

            for i, items in enumerate(classes.items()):
                if i == 0:
                    count_results_by_cluster_file.write('\t'.join([f"\t{items[0]}",
                                                        f"{items[1][0]}/{total_count_origin}",
                                                        f"{items[1][1]}/{total_count_inverse}" if reverse_matches else "-",
                                                        (f"{len(group_reps_ids[query])}") + 
                                                        (f"|{len(group_reps_ids[subject])}" if reverse_matches else "|-"),
                                                        (f"{len(group_alleles_ids[query])}" if reverse_matches else "-") +
                                                        (f"|{len(group_alleles_ids[subject])}"),
                                                        f"{frequency_in_genomes[query]}",
                                                        f"{frequency_in_genomes[subject]}\n"]))
                else:
                    count_results_by_cluster_file.write('\t'.join([f"\t\t{items[0]}",
                                                        f"{items[1][0]}/{total_count_origin}",
                                                        f"{items[1][1]}/{total_count_inverse}" if reverse_matches else "-",
                                                        "\n"]))
            count_results_by_cluster_file.write('\n')

def get_matches(all_relationships, clusters_to_keep, sorted_blast_dict):
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
    clusters_to_keep : dict
        A dictionary with classes as keys and lists of CDS or loci IDs to be kept as values.
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
    - The function first checks if `allele_ids` is provided to determine the mode of operation.
    - It uses utility functions like `itf.flatten_list` to simplify the structure of `all_relationships`
    and `itf.remove_by_regex` to clean up the IDs for matching.
    - The matching process accounts for whether entries are part of a joined cluster and adjusts the
    matching logic accordingly.
    - The function returns two dictionaries: one for general matches and one specifically for alleles,
    the latter being applicable only when `allele_ids` are provided.
    """
    is_matched = {}
    is_matched_alleles = {}

    relationships = itf.flatten_list(all_relationships.values())
    changed_ids = [[r[0], r[1].split('_')[0]] for r in relationships]
    had_matches = set([rep.split('_')[0] for rep in sorted_blast_dict])
    is_matched_alleles = {}
    for class_, entries in list(clusters_to_keep.items()):
        for entry in list(entries):
            if entry not in had_matches and not class_ == '1a':
                if entry == 'ERR5260641-protein1729':
                    print(entry)
                id_ = entry
                entry = [entry]
                is_matched.setdefault(id_, set([i[0] for i in changed_ids if i[1] in entry]))
                is_matched_alleles.setdefault(id_, set([i[1] 
                                                        for i in relationships 
                                                        if i[0] in is_matched[id_] 
                                                        and i[1].split('_')[0] in entry]))
    return is_matched, is_matched_alleles

def write_temp_loci(clusters_to_keep, not_included_cds, clusters, output_path):
    """
    This function wraps up the results for processing of the unclassified CDSs
    by writing FASTAs files for the possible new loci to include.
    
    Parameters
    ----------
    clusters_to_keep : dict
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
    drop_possible_loci : set
        Possible new loci that were dropped
    loci : dict
        Dict that contains the loci IDs and paths.
    groups_paths_old : dict
        The dictionary containing the old paths for the CDSs groups used 
        to cp instead of creating new FASTAs files.
    frequency_in_genomes : dict
        Dict that contains sum of frequency of that representatives cluster in the
        genomes of the schema.
    run_type : list
        What type of run to make.

    Returns
    -------
    groups_paths : dict
        Dict that contains as Key the ID of each group while the value is the
        path to the FASTA file that contains its nucleotide sequences.
    trans_dict_cds : dict
        Dict that contais the translations of all the CDSs inside the various
        groups.
    master_file_rep : str or None
        Path to the master file that contains all of the representative sequences.
    """

    def write_possible_new_loci(class_, cds_list, temp_fastas,
                                groups_paths, not_included_cds,
                                clusters):
        """
        Process each class and CDS list in clusters_to_keep.

        Parameters
        ----------
        class_ : str
            The class type.
        cds_list : list
            The list of CDSs.
        temp_fastas : str
            The path to the temp FASTAs folder.
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
            main_rep = cds
            cds_group_fasta_file = os.path.join(temp_fastas, main_rep + '.fasta')
            groups_paths[main_rep] = cds_group_fasta_file
            save_ids_index = {}
            if class_ != '1a':
                cds = [cds]
            else:
                cds = clusters_to_keep[class_][cds]
            index = 1
            # Write all of the alleles to the files.
            with open(cds_group_fasta_file, 'w') as fasta_file:
                for rep_id in cds:
                    cds_ids = clusters[rep_id]
                    for cds_id in cds_ids:
                        # Save the new ID to the dictionary where the old ID is the key.
                        save_ids_index[cds_id] = cds_id
                        # Write the allele to the file.
                        fasta_file.write(f">{index}\n{str(not_included_cds[cds_id])}\n")
                        index += 1

    temp_fastas_paths = {}
    print("Writing FASTA and additional files for possible new loci...")
    temp_fastas = os.path.join(output_path, 'temp_fastas')
    ff.create_directory(temp_fastas)
    # Process each class and CDS list in clusters_to_keep
    for class_, cds_list in clusters_to_keep.items():
        write_possible_new_loci(class_, cds_list, temp_fastas,
                                temp_fastas_paths, not_included_cds,
                                clusters)

    return temp_fastas_paths

def run_blasts(blast_db, cds_to_blast, reps_translation_dict,
               rep_paths_nuc, output_dir, constants, cpu, multi_fasta, run_type):
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
    multi_fasta : dict
       A dictionary used when the input FASTA files contain multiple CDSs, to ensure correct BLASTn
       execution.
    run_type : str
        A flag indicating what type of run to perform, can be cds_vs_cds, loci_vs_cds or loci_vs_loci.
        
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
    
    print("\nRunning BLASTn between cluster representatives vs cluster alleles..." if run_type == 'cds_vs_cds' else
          "\nRunning BLASTn between Schema representatives CDS clusters..." if run_type == 'loci_vs_cds' else
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
                res[1], constants[1], True, False, True, True, False)
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
        
        filename = itf.identify_string_in_dict_get_key(query_id, multi_fasta)
        if filename:
            blasts_to_run.setdefault(filename, set()).update(subjects_ids)
            seen_entries[filename] = set()
        else:
            filename = query_id
            seen_entries[filename] = set()
            blasts_to_run.setdefault(filename, set()).update(subjects_ids)
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
            
            _, self_score, _, _ = af.get_alignments_dict_from_blast_results(res[1], 0, True, True, True, True, False)
    
            # Save self-score
            self_score_dict[res[0]] = self_score
                            
            print(f"\rRunning BLASTp to calculate self-score for {res[0]: <{max_id_length}}", end='', flush=True)
            i += 1
    # Print newline
    print('\n')  
    
    print("Running BLASTp for representatives against cluster alleles..." if run_type == 'cds_vs_cds'
          else "Running BLASTp of schema representatives against cluster alleles..." if run_type == 'loci_vs_cds'
          else "Running BLASTp for schema representatives against schema alleles")
    # Run BLASTp between all BLASTn matches (rep vs all its BLASTn matches)  .      
    i = 1
    with concurrent.futures.ProcessPoolExecutor(max_workers=cpu) as executor:
        for res in executor.map(bf.run_blast_fastas_multiprocessing,
                                blastp_runs_to_do, 
                                repeat(get_blastp_exec),
                                repeat(blastp_results_folder),
                                repeat(rep_paths_prot),
                                rep_matches_prot.values()):
            
            filtered_alignments_dict, _, _, _ = af.get_alignments_dict_from_blast_results(res[1], 0, True, False, True, True, False)

            # Since BLAST may find several local aligments choose the largest one to calculate BSR.
            for query, subjects_dict in filtered_alignments_dict.items():
                for subject_id, results in subjects_dict.items():
                    #Highest score (First one)
                    subject_score = next(iter(results.values()))['score']
                    bsr_values[query].update({subject_id: bf.compute_bsr(subject_score, self_score_dict[res[0]])})
        
            print(f"\rRunning BLASTp for cluster representatives matches: {res[0]} - {i}/{total_blasts: <{max_id_length}}", end='', flush=True)
            i += 1

    return [representative_blast_results, representative_blast_results_coords_all,
            representative_blast_results_coords_pident, bsr_values, self_score_dict]

def write_processed_results_to_file(clusters_to_keep, representative_blast_results,
                                    classes_outcome, all_alleles, alleles, is_matched,
                                    is_matched_alleles, all_loci, output_path):
    """
    Write processed results to files in specified output directories.

    Parameters
    ----------
    clusters_to_keep : dict
        Dictionary containing clusters to keep, categorized by class.
    representative_blast_results : dict
        Dictionary containing representative BLAST results.
    classes_outcome : list
        List of class outcomes to process.
    all_alleles : dict
        Dictionary containing all alleles information.
    alleles : dict
        Dictionary containing specific alleles information.
    is_matched : dict
        Dictionary indicating if clusters are matched.
    is_matched_alleles : dict
        Dictionary containing matched alleles information.
    all_loci : bool
        Boolean indicating if all loci should be processed.
    output_path : str
        Path to the output directory where results will be written.

    Returns
    -------
    None
        This function does not return any value. It writes results to files.
    """
    # Create output directories
    blast_by_cluster_output = os.path.join(output_path, 'blast_by_cluster')
    ff.create_directory(blast_by_cluster_output)
    blast_results_by_class_output = os.path.join(output_path, 'blast_results_by_class')
    ff.create_directory(blast_results_by_class_output)

    add_group_column = True if 'loci_vs_cds' else False

    # Process clusters
    for class_, cds in clusters_to_keep.items():
        if class_ == 'Retained_not_matched_by_blastn':
            continue

        for cluster in cds if not isinstance(cds, dict) else cds.items():
            if isinstance(cds, dict):
                cluster_id = cluster[0]
                cluster = cluster[1]
                cluster_type = 'joined_cluster'
                cluster = cds[cluster_id]
            else:
                cluster_id = cluster
                cluster = [cluster]
                cluster_type = 'retained'
            
            is_cds = False
            if all_alleles:
                cluster_alleles = []
                for entry in cluster:
                    if alleles.get(entry):
                        cluster = alleles[entry]
                        cluster_type = 'CDS_cluster'
                        is_cds = True
                    else:
                        cluster_type = 'loci'
                        cluster_alleles += all_alleles[entry]
                if not is_cds:
                    cluster = cluster_alleles


            write_dict = {
                    query: {
                        subject: {id_: entry for id_, entry in entries.items()}
                        for subject, entries in subjects.items()
                    }
                    for query, subjects in representative_blast_results.items()
                    if query.split('_')[0] in cluster
                }
            
            if is_cds and class_ != '1a':
                queries = is_matched[cluster_id]
                cluster = is_matched_alleles[cluster_id]
                write_dict = {
                    query: {
                        subject: {id_: entry for id_, entry in entries.items()}
                        for subject, entries in subjects.items() if subject in cluster
                    }
                    for query, subjects in representative_blast_results.items()
                    if query in queries
                }

            report_file_path = os.path.join(blast_by_cluster_output, f"blast_{cluster_type}_{cluster_id}.tsv")
            alignment_dict_to_file(write_dict, report_file_path, 'w', add_group_column)

    # Process classes
    for class_ in classes_outcome:
        write_dict = {
            query: {
                subject: {id_: entry for id_, entry in entries.items() if entry['class'] == class_}
                for subject, entries in subjects.items()
            }
            for query, subjects in representative_blast_results.items()
        }
        report_file_path = os.path.join(blast_results_by_class_output, f"blastn_group_{class_}.tsv")
        alignment_dict_to_file(write_dict, report_file_path, 'w', add_group_column)

def extract_clusters_to_keep(classes_outcome, count_results_by_class, drop_mark):
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

    Returns
    -------
    clusters_to_keep : dict
        A dictionary with class identifiers as keys and lists of CDS identifiers or pairs of identifiers
        to be kept in each class.
    drop_possible_loci : set
        A set of CDS identifiers that are determined to be dropped based on their classification and
        presence in `drop_mark`.

    Notes
    -----
    - The function first initializes `clusters_to_keep` with empty lists for each class in `classes_outcome`.
    - It then iterates through `count_results_by_class` to assign CDS to the most appropriate class
    based on the provided outcomes.
    - Special handling is given to class '1a', where CDS pairs are clustered and indexed.
    - CDS marked in `drop_mark` and falling under certain classes are added to `drop_possible_loci` for exclusion.
    - The function uses utility functions like `itf.try_convert_to_type` for type conversion and
    `cf.cluster_by_ids` for clustering CDS pairs in class '1a'.
    """
    temp_keep = {}
    clusters_to_keep = {class_: [] for class_ in classes_outcome}
    drop_possible_loci = set()
    for ids, result in count_results_by_class.items():
        class_ = next(iter(result))
        [query, subject] = list(map(lambda x: itf.try_convert_to_type(x, int), ids.split('|')))
        if class_ == '1a':
            clusters_to_keep.setdefault('1a', []).append([query, subject])
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
            drop_possible_loci.add(itf.try_convert_to_type(keep, int))
        else:
            clusters_to_keep.setdefault(class_, []).append(keep)

    clusters_to_keep['1a'] = {i: list(values) for i, values in enumerate(cf.cluster_by_ids(clusters_to_keep['1a']), 1)}

    return clusters_to_keep, drop_possible_loci

def count_number_of_reps_and_alleles(clusters_to_keep, clusters, drop_possible_loci, group_reps_ids, group_alleles_ids):
    """
    Counts the number of representatives and alleles for each group in the given CDS clusters, excluding those in the drop set.

    Parameters
    ----------
    clusters_to_keep : dict
        Dictionary of CDS clusters to keep, organized by class and group.
    clusters : dict
        Dictionary mapping group IDs to their member CDS IDs.
    drop_possible_loci : set
        Set of group IDs to be excluded from the count.
    group_reps_ids : dict
        Dictionary to be updated with representative IDs for each group.
    group_alleles_ids : dict
        Dictionary to be updated with allele IDs for each group.

    Returns
    -------
    group_reps_ids : dict
        Dictionary where key is the CDS cluster ID and value is a set of representative IDs.
    group_alleles_ids : dict
        Dictionary where key is the CDS cluster ID and value is a set of allele IDs.
    """
    # Iterate over each class.
    for class_, cds_group in list(clusters_to_keep.items()):
        # Iterate over each group in class.
        for group in cds_group:
            if class_ == '1a':
                # Iterate over each representative in joined group.
                for cds in cds_group[group]:
                    if group_reps_ids.get(cds):
                        continue
                    group_reps_ids.setdefault(cds, set()).add(cds)
                    group_alleles_ids.setdefault(cds, set()).update(clusters[cds])
            elif group_reps_ids.get(group):
                continue
            else:
                group_reps_ids.setdefault(group, set()).add(group)
                group_alleles_ids.setdefault(group, set()).update(clusters[group])
    
    for id_ in drop_possible_loci:
        if id_ not in group_reps_ids:
            group_reps_ids.setdefault(id_, set()).add(id_)
            group_alleles_ids.setdefault(id_, set()).update(clusters[id_])

    return group_reps_ids, group_alleles_ids

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

def print_classifications_results(clusters_to_keep, drop_possible_loci, groups_paths_old, clusters, loci, run_type):
    """
    Prints the classification results based on the provided parameters.

    Parameters
    ----------
    clusters_to_keep : dict
        Dictionary containing CDS to keep, classified by their class type.
    drop_possible_loci : list
        List of possible loci dropped.
    groups_paths_old : dict
        The dictionary containing the old paths to the groups.
        Can be None.
    clusters : dict
        The dictionary containing the clusters.
    loci : bool
        If True, the analysis is based on loci.
    run_type : str
        The type of run, e.g., 'loci_vs_loci'.

    Returns
    -------
    cds_cases: dict
        Dictionary with CDS cases classified by their class type.
    loci_cases: dict
        Dictionary with loci cases classified by their class type.
    """
    def print_results(class_, count, printout, i):
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
                if run_type == 'loci_vs_loci':
                    print(f"\t\tOut of those {count} loci are recommended to be removed.")
                else:
                    print(f"\t\tOut of those {count} {'CDSs groups' if i== 0 else 'loci'}"
                        f" {'were removed from the analysis' if i == 0 else 'are recommended to be replaced with their matched CDS in the schema.'}")
            else:
                print(f"\t\tOut of those groups, {count} {'CDSs' if i == 0 else 'loci'} are classified as {class_} and were retained.")

    # If 'Retained_not_matched_by_blastn' exists in clusters_to_keep, remove it and store it separately
    Retained_not_matched_by_blastn = clusters_to_keep.pop('Retained_not_matched_by_blastn', None)

    # Display info about the results obtained from processing the classes.
    # Get the total number of CDS reps considered for classification.
    count_cases = {}
    loci_cases = {}
    cds_cases = {}
    if loci:
        # Iterate over classes and their associated CDS sets
        for class_, cds_set in clusters_to_keep.items():
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

        # Process drop_possible_loci in the same way as above
        loci_cases['dropped'] = [d for d in drop_possible_loci if d in loci]
        cds_cases['dropped'] = [d for d in drop_possible_loci if d not in loci]

    else:
        for class_, cds_set in clusters_to_keep.items():
            cds_cases[class_] = cds_set
            if class_ == '1a':
                count_cases[class_] = len(itf.flatten_list(cds_set.values()))
            else:
                count_cases[class_] = len(cds_set)
        cds_cases['dropped'] = drop_possible_loci
    # Check if loci is not empty
    if loci:
        for i, printout in enumerate([cds_cases, loci_cases]):
            if run_type == 'loci_vs_loci' and i == 0:
                continue
            total_loci = len(itf.flatten_list([i 
                                               for class_, i
                                               in printout.items()
                                               if class_ != '1a'])) + len(itf.flatten_list(clusters_to_keep['1a'].values()))

            print(f"Out of {len(groups_paths_old) if i==0 else len(loci)} {'CDSs groups' if i == 0 else 'loci'}:")
            print(f"\t{total_loci} {'CDSs' if i == 0 else 'loci'}"
                f" representatives had matches with BLASTn against the {'CDSs' if i == 1 else 'schema'}.")

            # Print the classification results
            for class_, group in printout.items():
                print_results(class_ ,len(group) if class_ != '1a' else len(itf.flatten_list(group.values())) ,printout, i)

            if i == 0:
                print(f"\t{len(groups_paths_old) - len(itf.flatten_list(printout.values()))}"
                    " didn't have any BLASTn matches so they were retained.\n")
    else:
        # Write info about the classification results.
        print(f"Out of {len(clusters)} clusters:")
        print(f"\t{sum(count_cases.values()) + len(drop_possible_loci)} CDS representatives had matches with BLASTn"
            f" which resulted in {len(itf.flatten_list(clusters_to_keep.values()))} groups")

        # Print the classification results
        for class_, count in count_cases.items():
            print_results(class_, count, clusters_to_keep, 0)

        print(f"\t\tOut of those {len(drop_possible_loci)} CDSs groups were removed from the analysis.")

        if Retained_not_matched_by_blastn:
            print(f"\t\t{len(Retained_not_matched_by_blastn)} didn't have any BLASTn matches so they were retained.")
            
            clusters_to_keep['Retained_not_matched_by_blastn'] = Retained_not_matched_by_blastn

    return cds_cases, loci_cases
