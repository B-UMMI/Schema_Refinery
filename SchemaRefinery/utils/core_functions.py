import os
import pandas as pd
import shutil
import time
from typing import Dict, Any, List, Tuple, Union, Optional, Set

try:
    from utils import (file_functions as ff,
                                        clustering_functions as cf,
                                        blast_functions as bf,
                                        alignments_functions as af,
                                        iterable_functions as itf,
                                        linux_functions as lf,
                                        graphical_functions as gf,
                                        pandas_functions as pf,
                                        sequence_functions as sf,
                                        Types as tp,
                                        print_functions as prf)
except ModuleNotFoundError:
    from SchemaRefinery.utils import (file_functions as ff,
                                        clustering_functions as cf,
                                        blast_functions as bf,
                                        alignments_functions as af,
                                        iterable_functions as itf,
                                        linux_functions as lf,
                                        graphical_functions as gf,
                                        pandas_functions as pf,
                                        sequence_functions as sf,
                                        Types as tp,
                                        print_functions as prf)

def alignment_dict_to_file(blast_results_dict: Dict[str, Dict[str, Dict[str, Dict[str, Any]]]], 
                           file_path: str, write_type: str) -> None:
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

    Returns
    -------
    None
        Writes the alignment data to the specified file based on the provided dictionary structure and parameters.

    Notes
    -----
    - It iterates over the nested dictionary structure to write each piece of alignment data to the file, formatting
    each row as tab-separated values.
    - This function is useful for exporting BLAST alignment results to a file for further analysis or reporting.
    """
    
    # Construct the header for the output file
    header: List[str] = ['Query\t',
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

    # Write or append to the file
    with open(file_path, write_type) as report_file:
        # Write the header only if the file is being created
        if write_type == 'w':
            report_file.write("".join(header))
        
        # Write all the alignment data
        for query_id, subjects in blast_results_dict.items():
            for subject_id, alignments in subjects.items():
                for alignment_id, alignment_data in alignments.items():
                    # Convert alignment data to a tab-separated string and write to the file
                    report_file.write('\t'.join(map(str, alignment_data.values())) + '\n')


# def sort_blast_results_by_classes(representative_blast_results: tp.BlastDict,
#                                   classes_outcome: Tuple[str, ...]
#                                   ) -> tp.BlastDict:
#     """
#     Sorts BLAST results by classes based on the alignment score.

#     This function organizes BLAST results into a sorted structure according to predefined classes.
#     It ensures that for each query, the results are grouped by the class of the alignment, prioritizing
#     the classes as specified in the `classes_outcome` list.

#     Parameters
#     ----------
#     representative_blast_results : Dict[str, Dict[str, List[Dict[str, Any]]]]
#         A dictionary where each key is a query identifier and each value is another dictionary.
#         The inner dictionary's keys are subject identifiers, and values are lists containing
#         details of the match, where the second element is a dictionary with the key 'class'
#         indicating the class of the alignment.
#     classes_outcome : Tuple[str, ...]
#         A list of possible classes outcomes to sort the BLAST results into. The order in this list
#         determines the priority of the classes when organizing the results.

#     Returns
#     -------
#     tp.BlastDict
#         A dictionary structured similarly to `representative_blast_results`, but sorted such that
#         all results for a given query are grouped by their class as determined by the highest
#         scoring alignment.

#     Notes
#     -----
#     - The function assumes that each match list in the values of `representative_blast_results`
#       contains at least one element, which is a dictionary with a 'class' key.
#     - It creates a temporary dictionary to first group results by class, then consolidates these
#       into the final sorted dictionary to be returned.
#     """
#     # Initialize the sorted dictionary
#     sorted_blast_dict: Dict[str, Dict[str, List[Dict[str, Any]]]] = {}

#     # Temporary dictionary to group results by class
#     temp_dict: Dict[str, Dict[str, Dict[str, List[Dict[str, Any]]]]] = {class_: {} for class_ in classes_outcome}

#     # Group results by class
#     for query, rep_blast_result in representative_blast_results.items():
#         for id_subject, matches in rep_blast_result.items():
#             # Get the class of the alignment with the highest score
#             class_ = matches[1]['classification']
#             if query not in temp_dict[class_]:
#                 temp_dict[class_][query] = {}
#             temp_dict[class_][query][id_subject] = matches

#     # Consolidate the grouped results into the final sorted dictionary
#     for class_ in classes_outcome:
#         for query, rep_blast_result in temp_dict[class_].items():
#             if query not in sorted_blast_dict:
#                 sorted_blast_dict[query] = {}
#             for id_subject, matches in rep_blast_result.items():
#                 sorted_blast_dict[query][id_subject] = matches

#     return sorted_blast_dict


def process_classes(representative_blast_results: tp.BlastDict,
                    representative_blastn_results: tp.BlastDict,
                    classes_outcome: Tuple[str, ...],
                    all_alleles: Dict[str, List[str]]) -> Tuple[
                    tp.ProcessedResults,
                    tp.CountResultsByClass,
                    tp.CountResultsByClassWithInverse,
                    tp.RepsAndAllelesIds,
                    List[str],
                    tp.AllRelationships,]:
    """
    Processes BLAST results to determine class-based relationships and counts.

    This function iterates through representative BLAST results to establish relationships
    between different coding sequences (CDS) and to count occurrences by class. It handles
    allele replacements, prioritizes classes based on a predefined order, and identifies
    important relationships between sequences.

    Parameters
    ----------
    representative_blast_results : tp.BlastDict
        A nested dictionary with the BLASTp results where the first key is the query sequence ID, the second key is
        the subject sequence ID, and the value is another dictionary containing match details
        including the class of the match.process_classes
    representative_blastn_results : tp.BlastDict
        A nested dictionary with the BLASTn results where the first key is the query sequence ID, the second key is
        the subject sequence ID, and the value is another dictionary containing match details
        including the class of the match.process_classes
    classes_outcome : Tuple[str, ...]
        A list of class identifiers ordered by priority. This order determines which classes are
        considered more significant when multiple matches for the same pair of sequences are found.
    all_alleles : Dict[str, List[str]]
        A dictionary mapping sequence IDs to their corresponding allele names. If provided, it is
        used to replace allele IDs with loci/CDS names in the processing.

    Returns
    -------
    processed_results : tp.ProcessedResults
        A dictionary containing processed results with keys formatted as "query|subject" and values being tuples
        containing information about the processed sequences, their class, relationships, and additional details.
    count_results_by_class : tp.CountResultsByClass
        A dictionary containing counts of results by class, with keys formatted as "query|subject" and values being
        dictionaries with class identifiers as keys and counts as values.
    count_results_by_class_with_inverse : tp.CountResultsByClassWithInverse
        A dictionary containing counts of results by class, including inverse matches, with keys formatted as
        "query|subject" and values being dictionaries with class identifiers as keys and counts as values.
    reps_and_alleles_ids : tp.RepsAndAllelesIds
        A dictionary mapping pairs of query and subject sequences to their unique loci/CDS IDs and alleles IDs.
    drop_mark : List[str]
        A list of sequences that were marked for dropping.
    all_relationships : tp.AllRelationships
        A dictionary containing all relationships between sequences, grouped by class.

    Notes
    -----
    - The function dynamically adjusts based on the presence of `all_alleles`, affecting how sequence IDs
    are replaced and processed.
    - It employs a complex logic to handle different scenarios based on class types and the presence or absence of
    alleles in the processed results, including handling allele replacements and determining the importance of
    relationships.
    """
    # Initialize variables
    count_results_by_class: tp.CountResultsByClass = {}
    count_results_by_class_with_inverse: tp.CountResultsByClassWithInverse = {}
    reps_and_alleles_ids: tp.RepsAndAllelesIds = {}
    processed_results: tp.ProcessedResults = {}
    all_relationships: tp.AllRelationships = {class_: [] for class_ in classes_outcome}
    drop_mark: List[str] = []
    inverse_match: List[str] = []
    files: List[str] = [representative_blast_results, representative_blastn_results]
    pattern: str = r'_(\d+)$'

    # Process the CDS to find what CDS to retain while also adding the relationships between different CDS with the BLAST results
    i = 0
    for file in files:
        i += 1
        for query, rep_blast_result in file.items():
            for id_subject, matches in rep_blast_result.items():
                class_: str = matches[1]['classification']
                all_relationships[class_].append([query, id_subject])
                ids_for_relationship: List[str] = [query, id_subject]
                query = itf.remove_by_regex(query, pattern)
                id_subject = itf.remove_by_regex(id_subject, pattern)
                new_query: str = query
                new_id_subject: str = id_subject

                strings: List[str] = [str(query), str(id_subject), class_]
                replaced_query: Optional[str] = itf.identify_string_in_dict_get_key(query, all_alleles)
                if replaced_query is not None:
                    new_query = replaced_query
                    strings[0] = new_query
                replaced_id_subject: Optional[str] = itf.identify_string_in_dict_get_key(id_subject, all_alleles)
                if replaced_id_subject is not None:
                    new_id_subject = replaced_id_subject
                    strings[1] = new_id_subject

                current_allele_class_index: int = classes_outcome.index(class_)
                # Check if the current loci were already processed
                if not processed_results.get(f"{new_query}|{new_id_subject}"):
                    run_next_step: bool = True
                # If those loci/CDS were already processed, check if the current class is better than the previous one
                elif current_allele_class_index < classes_outcome.index(processed_results[f"{new_query}|{new_id_subject}"][0]):
                    run_next_step = True
                # If not then skip the current alleles
                else:
                    run_next_step = False

                count_results_by_class.setdefault(f"{new_query}|{new_id_subject}", {})
                # Count the pairs in each class
                # Get the class of that pair
                if not count_results_by_class[f"{new_query}|{new_id_subject}"].get(class_):
                    # If the class has not yet be counted start the counter with 1
                    count_results_by_class[f"{new_query}|{new_id_subject}"].setdefault(class_, 1)
                else:
                    # Else add 1 to the counter
                    count_results_by_class[f"{new_query}|{new_id_subject}"][class_] += 1
                
                # Count the inversed pairs in each class
                if f"{new_query}|{new_id_subject}" not in inverse_match:
                    count_results_by_class_with_inverse.setdefault(f"{new_query}|{new_id_subject}", {})
                    inverse_match.append(f"{new_id_subject}|{new_query}")
                if f"{new_query}|{new_id_subject}" in inverse_match:
                    # If the class of the inversed pair has not yet been counted
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
                reps_and_alleles_ids.setdefault(f"{new_query}|{new_id_subject}", (set(), set()))
                if ids_for_relationship[0] not in reps_and_alleles_ids[f"{new_query}|{new_id_subject}"][0]:
                    reps_and_alleles_ids[f"{new_query}|{new_id_subject}"][0].add(ids_for_relationship[0])
                if ids_for_relationship[1] not in reps_and_alleles_ids[f"{new_query}|{new_id_subject}"][1]:
                    reps_and_alleles_ids[f"{new_query}|{new_id_subject}"][1].add(ids_for_relationship[1])
        
                if run_next_step:
                    # Set all None to run newly for this query/subject combination
                    processed_results[f"{new_query}|{new_id_subject}"] = (None, [], [], (new_query, new_id_subject), [])

                    if class_ in ['1b', '2a', '3a', '4a']:
                        blastn_entry: Dict[str, Any] = matches[list(matches.keys())[0]]
                        # Determine if the frequency of the query is greater than the subject.
                        is_frequency_greater: bool = blastn_entry['frequency_in_genomes_query_cds'] >= blastn_entry['frequency_in_genomes_subject_cds']
                        # Determine if the query or subject should be dropped.
                        query_or_subject: List[str] = [new_id_subject] if is_frequency_greater else [new_query]
                    else:
                        query_or_subject = []

                    # For the related_matches.tsv file.
                    if class_ not in ['4c','5']:
                        # Add asterisk to the query or subject that was dropped.
                        if class_ in ['1b', '2a', '3a', '4a']:
                            blastn_entry = matches[list(matches.keys())[0]]
                            # Determine if the frequency of the query is greater than the subject.
                            is_frequency_greater = blastn_entry['frequency_in_genomes_query_cds'] >= blastn_entry['frequency_in_genomes_subject_cds']
                            # Determine if the query or subject should be dropped.
                            dropped: str = new_id_subject if is_frequency_greater else new_query
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
        if i == 1:
            prf.print_message('BLASTp results processed.', 'info')
        if i == 2:
            prf.print_message('BLASTn results processed.', 'info')

    return processed_results, count_results_by_class, count_results_by_class_with_inverse, reps_and_alleles_ids, drop_mark, all_relationships


def extract_results(processed_results: tp.ProcessedResults, count_results_by_class: tp.CountResultsByClass, 
                    frequency_in_genomes: Dict[str, int], frequency_in_genomes_second_schema: Dict[str, int], merged_all_classes: tp.MergedAllClasses, 
                    dropped_loci_ids: Set[str], classes_outcome: Tuple[str]) -> Tuple[tp.RelatedClusters, tp.Recomendations, tp.Recomendations]:
    """
    Extracts and organizes results from process_classes.

    Parameters
    ----------
    processed_results : tp.ProcessedResults
        The processed results data.
    count_results_by_class : tp.CountResultsByClass
        A dictionary with counts of results by class.
    frequency_in_genomes : Dict[str, int]
        A dictionary containing the frequency of the query and subject in genomes.
    frequency_in_genomes_second_schema : Dict[str, int]
        A dictionary containing the frequency of the query and subject in genomes of the second schema.
    merged_all_classes : tp.MergedAllClasses
        A dictionary containing identifiers to check for a joined condition.
    dropped_loci_ids : List[str]
        A list of possible loci to drop.
    classes_outcome : Tuple[str]
        A list of class outcomes.

    Returns
    -------
    Tuple[Dict[str, List[Any]], Dict[str, Dict[str, Any]], Dict[str, Dict[str, Any]]]
        A tuple containing:
        - related_clusters: All relationships between loci and CDS.
        - recommendations: Dict that groups CDS/loci by ID and that contains strings to write in output file.
        - low_freq_recommendations: Dict that groups CDS/loci classified with 4c and 5 by ID and that contains strings to write in output file.
    
    Notes
    -----
    - The function iterates over `processed_results` to organize and cluster related CDS/loci based
      on their classification outcomes and the presence in specific clusters.
    - It uses helper functions like `cf.cluster_by_ids` for clustering and `itf.identify_string_in_dict_get_key`
      for identifying if a query or subject ID is present in the clusters.
    """
    def cluster_data(processed_results: tp.ProcessedResult) -> Dict[int, List[str]]:
        """
        Cluster data based on a specific key extracted from the processed results.

        Parameters
        ----------
        processed_results : Dict[str, Dict[str, Any]]
            A dictionary containing the processed results to be clustered.

        Returns
        -------
        Dict[int, List[str]]
            A dictionary where each key is an integer starting from 1, and each value is a cluster of data
            based on the extracted key, excluding entries with specific identifiers.
        """
        key_extractor = lambda v: v[3]
        condition = lambda v: v[0] not in ['4c', '5']

        return {i: cluster for i, cluster in enumerate(cf.cluster_by_ids([key_extractor(v) for v in processed_results.values() if condition(v)]), 1)}

    def choice_data(processed_results: tp.ProcessedResults, to_cluster_list: Dict[int, List[str]]) -> Dict[int, List[str]]:
        """
        Select and cluster data based on specific conditions and a list of identifiers to cluster.

        Parameters
        ----------
        processed_results : tp.ProcessedResults
            A dictionary containing the processed results to be selected and clustered.
        to_cluster_list : Dict[int, List[str]]
            A list of identifiers indicating which entries should be considered for clustering.

        Returns
        -------
        Dict[int, List[str]]
            A dictionary where each key is an integer starting from 1, and each value is a cluster of data
            selected based on the given conditions and the to_cluster_list.
        """
        # Extract the key (loci/CDS IDs) from the processed results
        key_extractor = lambda v: v[3]
        # Additional condition for the choice data (* means Dropped loci ID)
        additional_condition = lambda v: '*' in v[4][0] or '*' in v[4][1]

        return {i: cluster for i, cluster in enumerate(cf.cluster_by_ids([key_extractor(v) 
                                                                          for v in processed_results.values()
                                                                          if v[0] not in ['1a', '1b', '2a', '3a', '4a', '4c','5']
                                                                          or ((itf.identify_string_in_dict_get_key(v[3][0], to_cluster_list)
                                                                               and additional_condition(v))
                                                                               or (itf.identify_string_in_dict_get_key(v[3][1], to_cluster_list)
                                                                                   and additional_condition(v)))]), 1)}

        """
        Code to separate the choices into clusters based on the class of choice.

        classes_to_fetch_choice = ['1c', '2b', '3b', '4b', '6'] # Classes of choice
        choices_to_cluster = {keys: [] for keys in classes_to_fetch_choice}
        # For each choice class, select every entry ID that are part of the same to_cluster_list and share similiar class.
        for v in processed_results.values():
            for class_ in classes_to_fetch_choice:
                if v[0] == class_ or (itf.identify_string_in_dict_get_key(v[3][0], to_cluster_list) and additional_condition(v)) or (itf.identify_string_in_dict_get_key(v[3][1], to_cluster_list) and additional_condition(v)):
                    if not itf.identify_string_in_dict_get_key(v[3], choices_to_cluster) and not itf.identify_string_in_dict_get_key((v[3][1], v[3][0]), choices_to_cluster):
                        choices_to_cluster[class_].append(key_extractor(v))
                else:
                    continue

        clustered_choices = {}
        i = 0
        # Cluster the selected choices
        for cluster_choice in choices_to_cluster.values():
            clustered = cf.cluster_by_ids(cluster_choice)
            for cluster in clustered:
                clustered_choices[i] = cluster
                i += 1
                
        return clustered_choices
        """
    def process_id(id_: str, to_cluster_list: Dict[int, List[str]], merged_all_classes: tp.MergedAllClasses) -> Tuple[str, str, str]:
        """
        Process an identifier to check its presence in specific lists.

        Parameters
        ----------
        id_ : str
            The identifier to be processed.
        to_cluster_list : Dict[int, List[str]]
            A list of identifiers to check for the presence of id_.
        merged_all_classes : tp.MergedAllClasses
            A dictionary containing identifiers to check for a joined condition.

        Returns
        -------
        Tuple[str, bool, bool]
            A tuple containing the original id, a str indicating the key of the entry where the string is present, or None if not found,
            and a str indicating the key of the entry where the string is present in 1a.
        """
        present: str = itf.identify_string_in_dict_get_key(id_, to_cluster_list)
        joined_id: str = itf.identify_string_in_dict_get_key(id_, merged_all_classes['1a'])
        return id_, present, joined_id

    def check_in_recommendations(id_: str, joined_id: str, recommendations: tp.Recomendations, key: str, categories: List[str]) -> bool:
        """
        Check if an identifier or its joined form is present in a set of recommendations.

        Parameters
        ----------
        id_ : str
            The original identifier to check.
        joined_id : str
            The joined form of the identifier to check.
        recommendations : tp.Recomendations
            A dictionary of recommendations to search within.
        key : str
            The key within the recommendations to search under.
        categories : List[str]
            A list of categories to consider when searching in the recommendations.

        Returns
        -------
        bool
            True if the identifier or its joined form is present in the recommendations under the specified categories,
            False otherwise.
        """
        return any((joined_id or id_) in itf.flatten_list([v for k, v in recommendations[key].items() if cat in k]) for cat in categories)

    def add_to_recommendations(category: str, id_to_write: str, key: str, recommendations: tp.Recomendations, id_class: str, category_id: str = None) -> None:
        """
        Add an identifier to the recommendations under a specific category.

        Parameters
        ----------
        category : str
            The category under which to add the identifier.
        id_to_write : str
            The identifier to add to the recommendations.
        key : str
            The key within the recommendations to add under.
        recommendations : tp.Recomendations
            The recommendations dictionary to modify.
        joined_id : str, optional
            The joined form of the identifier, if applicable. Default is None.

        Returns
        -------
        None
            This function does not return any value but modifies the `recommendations` dictionary in place.
        """
        if category_id is not None:  # For joined or choice categories
            recommendations[key].setdefault(f'{category}_{category_id}', {})[id_to_write] = id_class
        else:  # For keep or drop categories
            recommendations[key].setdefault(category, {})[id_to_write] = id_class

    def order_dict_by_first_value(data: Dict[str, Tuple[str, List[str], List[str], Tuple[str, str], List[str]]], order: Tuple[str]) -> Dict[str, Tuple[str, List[str], List[str], Tuple[str, str], List[str]]]:
        """
        Order the dictionary by the first value of each tuple according to a specified order tuple.

        Parameters
        ----------
        data : Dict[str, Tuple[str, List[str], List[str], Tuple[str, str], List[str]]]
            The dictionary to be ordered. The dictionary values are tuples containing various elements.
        order : Tuple[str]
            The order tuple specifying the desired order of the first values in the dictionary tuples.

        Returns
        -------
        Dict[str, Tuple[str, List[str], List[str], Tuple[str, str], List[str]]]
            The ordered dictionary.
        """
        # Create a dictionary to map order values to their indices for sorting
        order_index = {value: index for index, value in enumerate(order)}
        
        # Sort the dictionary by the first value of each tuple according to the order tuple
        sorted_data = dict(sorted(data.items(), key=lambda item: order_index.get(item[1][0], len(order))))
        
        return sorted_data

    def find_both_values_in_dict_list(query_id: str, subject_id: str, choice: Dict[str, List[str]]) -> Union[str, None]:
        """
        Find the key where both query_id and subject_id are present in the same list.

        Parameters
        ----------
        query_id : str
            The query identifier to find.
        subject_id : str
            The subject identifier to find.
        choice : Dict[str, List[str]]
            The dictionary containing lists of identifiers.

        Returns
        -------
        str
            The key where both query_id and subject_id are present, or an empty string if not found.
        """
        for key, values in choice.items():
            if query_id in values and subject_id in values:
                return key
        return None

    # Initialize dictionaries to keep track of various results
    related_clusters: tp.RelatedClusters = {}  # To keep track of the related clusters
    recommendations: tp.Recomendations = {}  # To keep track of the recommendations
    dropped_match: Dict[str, List[List[str]]] = {}  # To keep track of the dropped matches

    # Filter the processed results by the order of classes
    processed_results = order_dict_by_first_value(processed_results, classes_outcome)

    # Normal run, where IDs are only loci or CDS original IDs.
    to_cluster_list: Dict[int, List[str]] = cluster_data(processed_results)  # All of the cluster of related CDS/loci by ID.
    choice: Dict[int, List[str]] = choice_data(processed_results, to_cluster_list)  # All of the possible cases of choice.

    for results in processed_results.values():
        
        if results[0] in ['4c', '5', 'Retained_not_matched_by_blastn']:
            continue
        # Get the IDs and check if they are present in the to_cluster_list
        query_id, query_present, joined_query_id = process_id(results[3][0], to_cluster_list, merged_all_classes)
        subject_id, subject_present, joined_subject_id = process_id(results[3][1], to_cluster_list, merged_all_classes)
        
        key: str = query_present if query_present else subject_present
        
        direct_match_info = f"{count_results_by_class[f'{results[3][0]}|{results[3][1]}'][results[0]]}/{sum(count_results_by_class[f'{results[3][0]}|{results[3][1]}'].values())}"
        if frequency_in_genomes_second_schema is not None:
            related_clusters.setdefault(key, []).append(results[4]
                                                        + [direct_match_info]
                                                        + [str(frequency_in_genomes[results[3][0]])]
                                                        + [str(frequency_in_genomes[results[3][1]])]
                                                        + [str(frequency_in_genomes_second_schema[results[3][0]])]
                                                        + [str(frequency_in_genomes_second_schema[results[3][1]])]
                                                        )
        else:
            related_clusters.setdefault(key, []).append(results[4]
                                                        + [direct_match_info]
                                                        + [str(frequency_in_genomes[results[3][0]])]
                                                        + [str(frequency_in_genomes[results[3][1]])]
                                                        )

        recommendations.setdefault(key, {})
        # Checks to make, so that we can add the ID to the right recommendations.
        if_same_joined: bool = (joined_query_id == joined_subject_id) if joined_query_id and joined_subject_id else False
        if_joined_query: bool = check_in_recommendations(query_id, joined_query_id, recommendations, key, ['Join'])
        if_joined_subject: bool = check_in_recommendations(subject_id, joined_subject_id, recommendations, key, ['Join'])
        if_query_in_choice: bool = check_in_recommendations(query_id, joined_query_id, recommendations, key, ['Choice'])
        if_subject_in_choice: bool = check_in_recommendations(subject_id, joined_subject_id, recommendations, key, ['Choice'])
        if_query_dropped: bool = (joined_query_id or query_id) in dropped_loci_ids
        if_subject_dropped: bool = (joined_subject_id or subject_id) in dropped_loci_ids

        choice_id: str = find_both_values_in_dict_list(query_id, subject_id, choice) # Find the choice ID (both query and subject in the same list)

        # What IDs to add to the Keep, Drop and Choice.
        query_to_write: str = joined_query_id or query_id
        subject_to_write: str = joined_subject_id or subject_id

        joined_query_to_write: str = query_id
        joined_subject_to_write: str = subject_id

        # Process the joined cases
        if results[0] == '1a':
            # If it is part of a joined cluster, add to the recommendations
            if joined_query_id is not None:
                add_to_recommendations('Join', joined_query_to_write, key, recommendations, results[0], joined_query_id)
            if joined_subject_id is not None:
                add_to_recommendations('Join', joined_subject_to_write, key, recommendations, results[0], joined_subject_id)
        # Process the choice cases
        elif results[0] in ['1c', '2b', '3b', '4b', '6']:
            # If it is not dropped and not the same joined cluster, add to the choice recommendations
            if not if_query_dropped and not if_subject_dropped and not if_same_joined:
                # If not already reported (other Joined elements already reported with Joined ID)
                if not find_both_values_in_dict_list(query_to_write, subject_to_write, recommendations[key]):
                    add_to_recommendations(f'Choice', query_to_write, key, recommendations, results[0], choice_id)
                    add_to_recommendations(f'Choice', subject_to_write, key, recommendations, results[0], choice_id)

        # Process cases where some ID is dropped
        elif results[0] in ['1b', '2a', '3a', '4a']:
            # If it is part of a joined cluster and it gets dropped, add to the recommendations
            if (joined_query_id and '*' in results[4][0]) or (joined_subject_id and '*' in results[4][1]) and not if_same_joined:
                # If not already reported (other Joined elements already reported with Joined ID)
                if not find_both_values_in_dict_list(query_to_write, subject_to_write, recommendations[key]):
                    add_to_recommendations(f'Choice', query_to_write, key, recommendations, results[0], choice_id)
                    add_to_recommendations(f'Choice', subject_to_write, key, recommendations, results[0], choice_id)

            # If it is not part of a joined cluster and it gets dropped, add to the recommendations
            if if_query_dropped:
                # If it is not part of a joined cluster and it gets dropped, add to the recommendations as Dropped
                if not if_joined_query and not if_query_in_choice:
                    add_to_recommendations('Drop', query_to_write, key, recommendations, results[0])
                    dropped_match.setdefault(key, []).append([query_to_write, subject_to_write, subject_to_write])
            # If it is not part of a joined cluster and it gets dropped, add to the recommendations
            elif if_subject_dropped:
                # If it is not part of a joined cluster and it gets dropped, add to the recommendations as Dropped
                if not if_joined_subject and not if_subject_in_choice:
                    add_to_recommendations('Drop', subject_to_write, key, recommendations, results[0])
                    dropped_match.setdefault(key, []).append([subject_to_write, query_to_write, query_to_write])

    # Add cases where some ID matched with dropped ID, we need to add the ID that matched with the ID that made the other match
    # to be Dropped. e.g x and y matched with x dropping and x also matched with z, then we need to make a choice between x and z.
    # Not being used currently
    """
    for key, matches in matched_with_dropped.items():
        if dropped_match.get(key):
            for dropped in dropped_match[key]:
                for match_ in matches:
                    if match_[2] in dropped:
                        add_to_recommendations('Choice', dropped[2], key, recommendations, match_[3])
    """

    sort_order: List[str] = ['Join', 'Choice', 'Drop', 'Add']
    recommendations = {k: {l[0]: l[1] for l in sorted(v.items(), key=lambda x: sort_order.index(x[0].split('_')[0]))} for k, v in recommendations.items()}

    low_freq_recommendations: tp.Recomendations = {}
    for results in processed_results.values():
        if results[0] in ['4c', '5']:
            query_id = results[3][0]
            subject_id = results[3][1]

            if query_id not in recommendations.keys():
                key: str = query_id
            elif subject_id not in recommendations.keys():
                key: str = subject_id
            low_freq_recommendations.setdefault(key, {})

            if query_id not in recommendations.keys():
                add_to_recommendations('Add', query_id, key, low_freq_recommendations, results[0])
            if subject_id not in recommendations.keys():
                add_to_recommendations('Add', subject_id, key, low_freq_recommendations, results[0])

    sort_order: List[str] = ['Join', 'Choice', 'Drop', 'Add']
    low_freq_recommendations = {k: {l[0]: l[1] for l in sorted(v.items(), key=lambda x: sort_order.index(x[0].split('_')[0]))} for k, v in low_freq_recommendations.items()}

    return related_clusters, recommendations, low_freq_recommendations


def write_recommendations_summary_results(to_blast_paths: Dict[str, str],
                                            related_clusters: tp.RelatedClusters,
                                            count_results_by_class_with_inverse: tp.CountResultsByClassWithInverse, 
                                            group_reps_ids: Dict[str, List[str]],
                                            group_alleles_ids: Dict[str, List[str]], 
                                            frequency_in_genomes: Dict[str, int],
                                            frequency_in_genomes_second_schema: Dict[str, int],
                                            recommendations: tp.Recomendations, 
                                            low_freq_recommendations: tp.Recomendations,
                                            reverse_matches: bool,
                                            classes_outcome: Tuple[str, ...],
                                            output_directory: str) -> Tuple[str, str, str, Dict[str, int]]:
    """
    Writes summary results of BLAST analysis to TSV files.

    This function generates two files: 'related_matches.tsv' and 'count_results_by_cluster.tsv'.
    The 'related_matches.tsv' file contains information about related clusters, while
    'count_results_by_cluster.tsv' details the count of results by cluster and class, including a total count.

    Parameters
    ----------
    related_clusters : tp.RelatedClusters
        A dictionary where each key is a cluster identifier and each value is a list of tuples.
        Each tuple represents a related match with its details.
    count_results_by_class_with_inverse : tp.CountResultsByClassWithInverse
        A dictionary where each key is a cluster identifier separated by '|', and each value is another dictionary.
        The inner dictionary's keys are class identifiers, and values are counts of results for that class.
    group_reps_ids : Dict[str, List[str]]
        A dictionary mapping sequence identifiers to their representative IDs.
    group_alleles_ids : Dict[str, List[str]]
        A dictionary mapping sequence identifiers to their allele IDs.
    frequency_in_genomes : Dict[str, int]
        A dictionary mapping sequence identifiers to their frequency in genomes.
    frequency_in_genomes_second_schema : Dict[str, int]
        A dictionary mapping sequence identifiers to their frequency in genomes in the second schema.
    recommendations : tp.Recommendations
        A dictionary containing recommendations for each cluster based on the classification of the results.
    low_freq_recommendations : tp.Recommendations
        A dictionary containing recommendations for each cluster based on the classification, of 4c and 5, of the results.
    reverse_matches : bool
        A flag indicating whether there are reverse matches.
    classes_outcome : Tuple[str, ...]
        A tuple of possible class outcomes to sort the BLAST results into.
    results_output : str
        The path to the directory where the output files will be saved.
    
    Returns
    -------
    related_matches_path : str
        Path to the output file with the matches relationships. 
    count_results_by_cluster_path : str
        Path to the count by clusters output file. 
    recommendations_file_path : str
        Path to the output file with the recommendations to every locus.
    count_classes_final : Dict[str, int]
        Dictionary with the number of loci in each class in the final recommendations.

    Notes
    -----
    - The 'related_matches.tsv' file is formatted such that each related match is written on a new line,
      with details separated by tabs. A blank line is added after each cluster's matches.
    - The 'count_results_by_cluster.tsv' file includes the cluster identifier, class identifier, count of results,
      and total count of results for the cluster, with each piece of information separated by tabs.
      A blank line is added after each cluster's information.
    """

    # Create recommendation file
    recommendations_file_path: str = os.path.join(output_directory, "recommendations.tsv")

    # Write loci groups by recommendation
    matched_loci: List[str] = []
    count_classes: Dict[str, int] = {} # To write the final stats report in the log file
    with open(recommendations_file_path, 'w') as recommendations_report_file:
        recommendations_report_file.write("Locus\tAction\tClass\n")
        for key, recommendation in recommendations.items():
            for category, loci in recommendation.items():
                if 'Drop' in category:
                    category=category
                else:
                    category = category.split('_', 1)[0]
                for locus, id_class in loci.items():
                    # Each gene can only have one action associated
                    if locus not in matched_loci:
                        recommendations_report_file.write(f"{locus}\t{category}\t{id_class}\n")
                        matched_loci.append(locus)
                        if id_class not in count_classes.keys():
                            count_classes.update({id_class: 1})
                        else:
                            count_classes[id_class] += 1
                    else:
                        continue
            recommendations_report_file.write("#\n")
        for key, recommendation in low_freq_recommendations.items():
            for category, loci in recommendation.items():
                for locus, id_class in loci.items():
                    # Each gene can only have one action associated
                    if locus not in matched_loci:
                        recommendations_report_file.write(f"{locus}\t{category}\t{id_class}\n")
                        matched_loci.append(locus)
                        if id_class not in count_classes.keys():
                            count_classes.update({id_class: 1})
                        else:
                            count_classes[id_class] += 1
                    else:
                        continue
        # Add the loci that had no action needed to the output file with the action 'Add' and class '7'
        for loci, loci_path in to_blast_paths.items():
            loci_id = ff.file_basename(loci).split('.')[0]
            if loci_id not in matched_loci:
               recommendations_report_file.write(f"{loci_id}\tAdd\t7\n") 
    
    order_index = {value: index for index, value in enumerate(classes_outcome)}
    count_classes_final = dict(sorted(count_classes.items(), key=lambda item: order_index.get(item[0], len(classes_outcome))))
            
    # Add the reverse matches to the related clusters
    reported_cases: Dict[str, List[Tuple[str, str]]] = {}
    for key, related in list(related_clusters.items()):
        for index, r in enumerate(list(related)):
            if reverse_matches:
                r.insert(4, '-')
                r.insert(5, '-')
            [query, subject] = [itf.remove_by_regex(i, r"\*") for i in r[:2]]
            
            # Add the total count of the representatives and alleles
            representatives_count_query: str = f"{len(group_reps_ids[query])}"
            representatives_count_subject: str = f"{len(group_reps_ids[subject])}" if reverse_matches else "|-"
            allele_count_query: str = f"{len(group_alleles_ids[query])}" if reverse_matches else "-"
            allele_count_subject: str = f"{len(group_alleles_ids[subject])}"
            representatives_count: str = f"{representatives_count_query}|{representatives_count_subject}"
            allele_count: str = f"{allele_count_query}|{allele_count_subject}"
            
            r.extend([representatives_count, allele_count])
            if (query, subject) not in itf.flatten_list(reported_cases.values()):
                reported_cases.setdefault(key, []).append((subject, query))
            elif reverse_matches:
                sublist_index: int = itf.find_sublist_index([[itf.remove_by_regex(i, r"\*")
                                                              for i in l[:2]]
                                                              for l in related_clusters[key]], [subject, query])
                insert: str = r[2] if r[2] is not None else '-'
                related[sublist_index][4] = insert
                insert = r[3] if r[3] is not None else '-'
                related[sublist_index][5] = insert
                related.remove(r)

    # Write the results to the output files
    related_matches_path: str = os.path.join(output_directory, "related_matches.tsv")
    with open(related_matches_path, 'w') as related_matches_file:
        related_matches_file.write("Query\tSubject\tClass\tClass_count" +
                                    ("\tInverse_class\tInverse_class_count" if reverse_matches else "") +
                                    "\tFrequency_in_genomes_query\tFrequency_in_genomes_subject"
                                    "\tAlleles_used_to_blast_count\tAlleles_blasted_against_count\n")
        for related in related_clusters.values():
            for r in related:
                related_matches_file.write('\t'.join(str(item) for item in r) + '\n')

            related_matches_file.write('#\n')
    
    # Write the count results to the output file
    tab: str = '\t'
    count_results_by_cluster_path: str = os.path.join(output_directory, "count_results_by_cluster.tsv")
    with open(count_results_by_cluster_path, 'w') as count_results_by_cluster_file:
        # for run_mode == schema_vs_schema
        if frequency_in_genomes_second_schema is not None:
            count_results_by_cluster_file.write("Query"
                                                "\tSubject"
                                                f"\t{tab.join(classes_outcome)}"
                                                "\tAlleles_used_to_blast_count"
                                                "\tAlleles_blasted_against_count"
                                                "\tFrequency_in_genomes_query"
                                                "\tFrequency_in_genomes_subject"
                                                "\tFrequency_in_genomes_second_schema_query"
                                                "\tFrequency_in_genomes_second_schema_subject\n")
        else:
            count_results_by_cluster_file.write("Query"
                                                "\tSubject"
                                                f"\t{tab.join(classes_outcome)}"
                                                "\tAlleles_used_to_blast_count"
                                                "\tAlleles_blasted_against_count"
                                                "\tFrequency_in_genomes_query"
                                                "\tFrequency_in_genomes_subject\n")
        for id_, classes in count_results_by_class_with_inverse.items():
            total_count_origin: int = sum([i[0] for i in classes.values() if i[0] != '-'])
            total_count_inverse: int = sum([i[1] for i in classes.values() if i[1] != '-'])
            query = itf.try_convert_to_type(id_.split('|')[0], int)
            subject = itf.try_convert_to_type(id_.split('|')[1], int)
            
            representatives_count_query = f"{len(group_reps_ids[query])}"
            representatives_count_subject = f"{len(group_reps_ids[subject])}" if reverse_matches else "|-"
            allele_count_query = f"{len(group_alleles_ids[query])}" if reverse_matches else "-"
            allele_count_subject = f"{len(group_alleles_ids[subject])}"
            representatives_count = f"{representatives_count_query}|{representatives_count_subject}"
            allele_count = f"{allele_count_query}|{allele_count_subject}"
            
            query_frequency: int = frequency_in_genomes[query]
            subject_frequency: int = frequency_in_genomes[subject]
            # for run_mode == schema_vs_schema
            if frequency_in_genomes_second_schema is not None:
                query_frequency_ss: int = frequency_in_genomes_second_schema[query]
                subject_frequency_ss: int = frequency_in_genomes_second_schema[subject]

            # Start constructing the line
            line_parts = [query, subject]

            for class_outcome in classes_outcome:
                class_value = classes.get(class_outcome, '-')
                if class_value == '-':
                    line_parts.append('-')
                else:
                    query_class_count = class_value[0] if class_value[0] != '-' else '-'
                    subject_class_count = class_value[1] if class_value[1] != '-' else '-'
                    line_parts.append(f"{query_class_count}|{total_count_origin}|{subject_class_count}|{total_count_inverse}")

            # Add summary counts and frequencies
            line_parts.append(representatives_count)
            line_parts.append(allele_count)
            line_parts.append(str(query_frequency))
            line_parts.append(str(subject_frequency))

            if frequency_in_genomes_second_schema is not None:
                line_parts.append(str(query_frequency_ss))
                line_parts.append(str(subject_frequency_ss))

            count_results_by_cluster_file.write('\t'.join(line_parts) + '\n')

    return (related_matches_path, count_results_by_cluster_path, recommendations_file_path, count_classes_final)


def get_matches(all_relationships: tp.AllRelationships, merged_all_classes: tp.MergedAllClasses, 
                sorted_blast_dict: tp.BlastDict) -> Tuple[Dict[str, Set[str]], Union[Dict[str, Set[str]], None]]:
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
    all_relationships : tp.AllRelationships
        A dictionary containing all relationships between loci and alleles or CDS, with loci as keys
        and lists of related alleles or CDS as values.
    merged_all_classes : tp.MergedAllClasses
        A dictionary with classes as keys and lists of CDS or loci IDs to be kept as values.
    sorted_blast_dict : Dict[str, Any]
        A dictionary containing sorted BLAST results, used to identify loci that have matches.

    Returns
    -------
    is_matched : Dict[str, Set[str]]
        A dictionary with loci or CDS IDs as keys and sets of matched loci IDs as values, indicating
        successful matches.
    is_matched_alleles : Dict[str, Set[str]] or None
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
    # Initialize dictionaries to store matches
    is_matched: Dict[str, Set[str]] = {}
    is_matched_alleles: Dict[str, Set[str]] = {}

    # Flatten the list of all relationships
    relationships: List[Tuple[str, str]] = itf.flatten_list(all_relationships.values())
    
    # Change IDs to remove suffixes
    changed_ids: List[List[str]] = [[r[0], r[1].split('_')[0]] for r in relationships]
    
    # Identify loci that had matches in the sorted BLAST results
    had_matches: Set[str] = set([rep.split('_')[0] for rep in sorted_blast_dict])
    
    # Iterate over clusters to keep and determine matches
    for class_, entries in list(merged_all_classes.items()):
        for entry in list(entries):
            if entry not in had_matches and class_ != '1a':
                id_: str = entry
                entry_list: List[str] = [entry]
                
                # Determine general matches
                is_matched.setdefault(id_, set([i[0] for i in changed_ids if i[1] in entry_list]))
                
                # Determine allele-specific matches
                is_matched_alleles.setdefault(id_, set([i[1] 
                                                        for i in relationships 
                                                        if i[0] in is_matched[id_] 
                                                        and i[1].split('_')[0] in entry_list]))
    
    return is_matched, is_matched_alleles


def write_processed_results_to_file(merged_all_classes: tp.MergedAllClasses, 
                                    representative_blast_results: tp.BlastDict,
                                    classes_outcome: List[str], all_alleles: Dict[str, List[str]], 
                                    is_matched: Dict[str, Set[str]], is_matched_alleles: Dict[str, Set[str]], 
                                    output_path: str) -> None:
    """
    Write processed results to files in specified output directories.

    Parameters
    ----------
    merged_all_classes : tp.MergedAllClasses
        Dictionary containing clusters to keep, categorized by class.
    representative_blast_results : tp.BlastDict
        Dictionary containing representative BLAST results.
    classes_outcome : List[str]
        List of class outcomes to process.
    all_alleles : Dict[str, List[str]]
        Dictionary containing all alleles information.
    is_matched : Dict[str, Set[str]]
        Dictionary indicating if clusters are matched.
    is_matched_alleles : Dict[str, Set[str]]
        Dictionary containing matched alleles information.
    output_path : str
        Path to the output directory where results will be written.

    Returns
    -------
    None
        This function does not return any value. It writes results to files.
    """
    # Create output directories
    blast_by_cluster_output: str = os.path.join(output_path, 'blast_by_cluster')
    ff.create_directory(blast_by_cluster_output)
    blast_results_by_class_output: str = os.path.join(output_path, 'blast_results_by_class')
    ff.create_directory(blast_results_by_class_output)

    # Process clusters
    for class_, cds in merged_all_classes.items():
        # Skip clusters that are not matched
        if class_ == 'Retained_not_matched_by_blastn':
            continue
        # Process clusters
        for cluster in cds if not isinstance(cds, dict) else cds.items():
            # Determine cluster type
            if isinstance(cds, dict):
                cluster_id: str = cluster[0]
                cluster = cluster[1]
                cluster_type: str = 'joined_cluster'
                cluster = cds[cluster_id]
            else:
                cluster_id = cluster
                cluster = [cluster]
                cluster_type = 'retained'
            # Get all alleles in the cluster
            cluster_alleles: List[str] = []
            for entry in cluster:
                cluster_alleles += all_alleles[entry]
            # Write results
            write_dict: Dict[str, Dict[str, Dict[str, Any]]] = {
                query: {
                    subject: {id_: entry for id_, entry in entries.items()}
                    for subject, entries in subjects.items()
                }
                for query, subjects in representative_blast_results.items()
                if query.split('_')[0] in cluster
            }
            # If the cluster is matched and not joined, write the results
            if is_matched.get(cluster_id) and class_ != '1a':
                queries: Set[str] = is_matched[cluster_id]
                cluster = is_matched_alleles[cluster_id]
                write_dict = {
                    query: {
                        subject: {id_: entry for id_, entry in entries.items()}
                        for subject, entries in subjects.items() if subject in cluster
                    }
                    for query, subjects in representative_blast_results.items()
                    if query in queries
                }

            report_file_path: str = os.path.join(blast_by_cluster_output, f"blast_{cluster_type}_{cluster_id}.tsv")
            alignment_dict_to_file(write_dict, report_file_path, 'w')

    # Process classes
    for class_ in classes_outcome:
        write_dict = {
            query: {
                subject: {id_: entry for id_, entry in entries.items() if entry['classification'] == class_}
                for subject, entries in subjects.items()
            }
            for query, subjects in representative_blast_results.items()
        }
        report_file_path = os.path.join(blast_results_by_class_output, f"blastn_group_{class_}.tsv")
        alignment_dict_to_file(write_dict, report_file_path, 'w')


def extract_clusters_to_keep(classes_outcome: Tuple[str], count_results_by_class: tp.CountResultsByClass, 
                            drop_mark: List[str]
                            ) -> Tuple[Dict[int, List[str]], tp.MergedAllClasses, Set[str]]:
    """
    Extracts and organizes CDS (Coding Sequences) to keep based on classification outcomes.

    This function processes BLAST results to determine which coding sequences (CDS) should
    be retained for further analysis based on their classification outcomes. It organizes
    CDS into categories, prioritizes them according to a predefined order of classes, and
    identifies sequences to be dropped.

    Parameters
    ----------
    classes_outcome : List[str]
        An ordered list of class identifiers that determine the priority of classes for keeping CDS.
    count_results_by_class : tp.CountResultsByClass
        A dictionary where keys are concatenated query and subject IDs separated by '|', and values
        are dictionaries with class identifiers as keys and counts as values.
    drop_mark : Set[str]
        A set of identifiers that are marked for dropping based on previous criteria.

    Returns
    -------
    clusters_to_keep_1a: Dict[int, List[str]]
        A dictionary with class 1a cases.
    clusters_to_keep : tp.MergedAllClasses
        A dictionary with class identifiers as keys and lists of CDS identifiers or pairs of identifiers
        to be kept in each class.
    drop_possible_loci : Set[str]
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
    # Temporary dictionary to keep track of the best class for each CDS
    temp_keep: Dict[str, str] = {}
    
    # Initialize clusters_to_keep with empty lists for each class in classes_outcome
    clusters_to_keep: tp.MergedAllClasses = {class_: [] for class_ in classes_outcome if class_ != '1a'}
    temp_clusters_1a: List[List[str]] = []
    # Initialize a set to keep track of CDS to be dropped
    drop_possible_loci: Set[str] = set()
    
    # Iterate through count_results_by_class to assign CDS to the most appropriate class
    for ids, result in count_results_by_class.items():
        class_: str = next(iter(result))
        query: str
        subject: str
        query, subject = ids.split('|')
        
        # Special handling for class '1a'
        if class_ == '1a':
            temp_clusters_1a.append([query, subject])
        
        # Update temp_keep with the best class for each CDS
        if not temp_keep.get(query):
            temp_keep[query] = class_
        elif classes_outcome.index(class_) < classes_outcome.index(temp_keep[query]):
            temp_keep[query] = class_
        
        if not temp_keep.get(subject):
            temp_keep[subject] = class_
        elif classes_outcome.index(class_) < classes_outcome.index(temp_keep[subject]):
            temp_keep[subject] = class_
    
    # Iterate through temp_keep to finalize clusters_to_keep and drop_possible_loci
    for keep, class_ in temp_keep.items():
        if class_ == '1a':
            continue
        if keep in drop_mark and class_ in ['1b', '2a', '3a', '4a']:
            drop_possible_loci.add(itf.try_convert_to_type(keep, int))
        else:
            clusters_to_keep.setdefault(class_, []).append(keep)
    
    # Cluster and index CDS pairs in class '1a'
    clusters_to_keep_1a: Dict[int, List[str]] = {i: list(values) for i, values in enumerate(cf.cluster_by_ids(temp_clusters_1a), 1)}
    
    return clusters_to_keep_1a, clusters_to_keep, drop_possible_loci


def count_number_of_reps_and_alleles(merged_all_classes: tp.MergedAllClasses, 
                                    clusters: Dict[str, List[str]], 
                                    drop_possible_loci: Set[str], group_reps_ids: Dict[str, List[str]], 
                                    group_alleles_ids: Dict[str, List[str]]) -> Tuple[Dict[str, List[str]], Dict[str, List[str]]]:
    """
    Counts the number of representatives and alleles for each group in the given CDS clusters, excluding those in
    the drop set.

    This function iterates through the clusters to keep and counts the number of representatives and alleles
    for each group, updating the provided dictionaries with the results. It excludes groups that are marked
    for dropping.

    Parameters
    ----------
    clusters_to_keep : Dict[str, Dict[int, List[int]]]
        Dictionary of CDS clusters to keep, organized by class and group.
    clusters : Dict[int, List[int]]
        Dictionary mapping group IDs to their member CDS IDs.
    drop_possible_loci : Set[int]
        Set of group IDs to be excluded from the count.
    group_reps_ids : Dict[int, Set[int]]
        Dictionary to be updated with representative IDs for each group.
    group_alleles_ids : Dict[int, Set[int]]
        Dictionary to be updated with allele IDs for each group.

    Returns
    -------
    Tuple[Dict[int, Set[int]], Dict[int, Set[int]]]
        Updated dictionaries where the key is the CDS cluster ID and the value is a set of representative IDs
        and allele IDs, respectively.

    Notes
    -----
    - The function first iterates through the clusters to keep, updating the representative and allele IDs
      for each group.
    - It then iterates through the drop set to ensure that any groups marked for dropping are also included
      in the dictionaries.
    """
    def add_ids(cds: str):
        """
        Helper function to add representative and allele IDs to the respective dictionaries.
        """
        group_reps_ids.setdefault(cds, set()).add(cds)
        group_alleles_ids.setdefault(cds, set()).add(cds)

    # Iterate over each class in merged_all_classes
    for class_, cds_group in merged_all_classes.items():
        # Iterate over each group in the class
        for group in cds_group:
            if class_ == '1a':
                # Special handling for class '1a'
                for cds in cds_group[group]:
                    if group_reps_ids.get(cds):
                        # Skip if the CDS is already in group_reps_ids
                        continue
                    add_ids(cds)
            elif group_reps_ids.get(group):
                # Skip if the group is already in group_reps_ids
                continue
            else:
                add_ids(group)
    
    # Iterate over the drop set to ensure groups marked for dropping are included
    for id_ in drop_possible_loci:
        if id_ not in group_reps_ids:
            # Add the group ID to group_reps_ids and group_alleles_ids if not already present
            group_reps_ids.setdefault(id_, set()).add(id_)
            group_alleles_ids.setdefault(id_, set()).update(clusters[id_])

    return group_reps_ids, group_alleles_ids


def print_classifications_results(merged_all_classes: tp.MergedAllClasses, drop_possible_loci: List[int], 
                                  to_blast_paths: Dict[str, str], clusters: Dict[str, Any], count_classes_final: Dict[str, int]) -> None:
    """
    Prints the classification results based on the provided parameters.

    Parameters
    ----------
    merged_all_classes : tp.MergedAllClasses
        Dictionary containing CDS to keep, classified by their class type.
    drop_possible_loci : List[int]
        List of possible loci dropped.
    to_blast_paths : Dict[str, str]
        Path to BLAST.
    clusters : Dict[str, Any]
        The dictionary containing the clusters.
    moved_recs : Dict[str, List[Set[str]]]
        The dictionary of which loci in each class got the recommendation they were expected to or were moved, for exemple 'Drop' instead of 'Choice'.

    Returns
    -------
    None
        Prints the classification results to stdout.

    Notes
    -----
    - The function first processes the `clusters_to_keep` dictionary to count the number of CDS
      representatives for each class.
    - It then prints the results, including the number of groups classified under each class and
      any recommendations for verification.
    - If there are any retained groups not matched by BLASTn, it handles them separately.
    """
    def print_results(class_: str, count: int, printout: Dict[str, Any]) -> None:
        """
        Prints the classification results based on the class type.

        Parameters
        ----------
        class_ : str
            The class type.
        count : int
            The count of groups.
        printout : Dict[str, Any]
            The dictionary containing printout information.
        moved_recs : Dict[str, List[Set[str]]]
            The dictionary of which loci in each class got the recommendation they were expected to or were moved, for exemple 'Drop' instead of 'Choice'.

        Returns
        -------
        None
            Prints the classification results to stdout.

        Notes
        -----
        - The function prints different messages based on the class type.
        - It provides recommendations for verification for certain classes.
        """
        if count > 0:
            if class_ in ['2b', '4b']:
                prf.print_message(f"\t\t{count} loci are classified as {class_} and were retained"
                            " but it is recommended to verify them as they may be contained or contain partially inside"
                            " their BLAST match.", None)
            elif class_ == '1a':
                prf.print_message(f"\t\t{count} loci are classified as {class_}"
                            f" and are contained in {len(printout['1a'])} joined groups that were retained.", None)
            elif class_ in ['4c', '5']:
                prf.print_message(f"\t\t{count} loci are classified as {class_}. These will be added to recommendations with 'Add'.", None)
            else:
                prf.print_message(f"\t\t{count} loci are classified as {class_}.", None)

    # If 'Retained_not_matched_by_blastn' exists in clusters_to_keep, remove it and store it separately
    retained_not_matched_by_blastn: Optional[Any] = merged_all_classes.pop('Retained_not_matched_by_blastn', None)

    # Display info about the results obtained from processing the classes.
    # Get the total number of CDS reps considered for classification.
    count_cases: Dict[str, int] = {}
    # Check if loci is not empty
    total_loci: int = sum(count_cases.values()) + sum(count_classes_final.values())
    prf.print_message(f"Out of {len(to_blast_paths)}:", None)
    prf.print_message(f"\t{total_loci} representatives had matches with BLASTn against the schema DNA sequences.", None)
    for class_, count in count_classes_final.items():
        print_results(class_, count, merged_all_classes)
    for class_, count in count_cases.items():
        print_results(class_, count, merged_all_classes)
    prf.print_message(f"\tOut of those {len(to_blast_paths.values()) - sum(count_cases.values()) - sum(count_classes_final.values())} didn't have any matches", None)

    if retained_not_matched_by_blastn:
        merged_all_classes['Retained_not_matched_by_blastn'] = retained_not_matched_by_blastn


def add_cds_to_dropped_cds(drop_possible_loci: Set[str], dropped_cds: Dict[str, str], clusters_to_keep: Dict[str, Any],
                           clusters_to_keep_1a: Dict[str, List[str]],
                            clusters: Dict[str, List[str]], reason: str, processed_drop: List[str]) -> None:
    """
    Adds CDS to the dropped CDS list based on the provided parameters.

    This function processes a list of possible loci to be dropped and updates the dropped CDS list
    with the provided reason. It ensures that each CDS is processed only once and handles both
    individual CDS and joined groups of CDS.

    Parameters
    ----------
    drop_possible_loci : List[int]
        List of possible loci dropped.
    dropped_cds : Dict[int, str]
        Dictionary to store dropped CDS with their reasons.
    clusters_to_keep : Dict[str, Any]
        Dictionary containing CDS to keep, classified by their class type.
    clusters : Dict[int, List[int]]
        Dictionary containing clusters of CDS.
    reason : str
        Reason for dropping the CDS.
    processed_drop : List[int]
        List of already processed drop IDs.

    Returns
    -------
    None
        The function updates the `dropped_cds` dictionary in place.

    Notes
    -----
    - The function first checks if each drop ID has already been processed to avoid duplication.
    - It identifies whether the drop ID belongs to a joined group (class '1a') or another class.
    - If the drop ID belongs to a joined group, it processes all elements of the group.
    - It updates the `dropped_cds` dictionary with the reason for dropping each CDS.
    """
    # Iterate over a copy of the drop_possible_loci list to avoid modifying the list during iteration
    for drop_id in drop_possible_loci.copy():
        # Skip if the drop ID has already been processed
        if drop_id in processed_drop:
            continue
        else:
            processed_drop.append(drop_id)

        # Check if the drop ID belongs to a joined group (class '1a')
        class_1a_id: Optional[str] = itf.identify_string_in_dict_get_key(drop_id, clusters_to_keep_1a)
        if class_1a_id:
            # Add all elements of the joined group to the dropped list
            process_ids: List[str] = clusters_to_keep_1a[class_1a_id]
            del clusters_to_keep_1a[drop_id]
        else:
            # Check if the drop ID belongs to another class
            class_: Optional[str] = itf.identify_string_in_dict_get_key(drop_id, {key: value for key, value in clusters_to_keep.items() if key != '1a'})
            if class_:
                clusters_to_keep[class_].remove(drop_id)

        # Update the dropped_cds dictionary with the reason for dropping each CDS
        if class_1a_id:
            for rep_id in process_ids:
                for cds_id in clusters[rep_id]:
                    dropped_cds[cds_id] = reason
        else:
            for cds_id in clusters[drop_id]:
                dropped_cds[cds_id] = reason


def write_dropped_possible_new_loci_to_file(drop_possible_loci: Set[str], dropped_cds: Dict[str, str], 
                                            output_directory: str) -> str:
    """
    Write the dropped possible new loci to a file with the reasons for dropping them.

    This function writes the IDs of possible new loci that should be dropped, along with the reasons
    for dropping them, to a TSV file in the specified output directory.

    Parameters
    ----------
    drop_possible_loci : Set[str]
        A set of possible new loci IDs that should be dropped.
    dropped_cds : Dict[str, str]
        A dictionary where keys are CDS (Coding Sequences) IDs and values are the reasons for dropping them.
    results_output : str
        The path to the directory where the output file will be saved.

    Returns
    -------
    None
        The function writes the results to a file and does not return any value.

    Notes
    -----
    - The function first constructs the output file path.
    - It then creates a dictionary mapping locus IDs to their drop reasons by extracting the locus ID
      from the CDS ID.
    - Finally, it writes the locus IDs and their drop reasons to the output file.
    """
    # Construct the output file path
    drop_possible_loci_output: str = os.path.join(output_directory, 'drop_loci_reason.tsv')
    
    # Create a dictionary mapping locus IDs to their drop reasons
    locus_drop_reason: Dict[str, str] = {cds.rsplit('_', 1)[0]: reason 
                                         for cds, reason in dropped_cds.items() if '_' in cds}
    
    # Write the locus IDs and their drop reasons to the output file
    with open(drop_possible_loci_output, 'w') as drop_possible_loci_file:
        drop_possible_loci_file.write('Possible_new_loci_ID\tDrop_Reason\n')
        for locus in drop_possible_loci:
            drop_possible_loci_file.write(f"{locus}\t{locus_drop_reason[locus]}\n")
    
    return drop_possible_loci_output


def prepare_loci(schema_folder: str,
                 constants: List[Any],
                 results_output: str) -> Tuple[
                     Dict[str, str], 
                     str, 
                     Dict[str, str],
                     Dict[str, str],  
                     Dict[str, List[str]],  
                     Dict[str, List[str]], 
                     Dict[str, List[str]],
                     Dict[str, str],
                     Dict[str, int],
                     Dict[str, int]]:
    """
    Process new loci by translating sequences, counting frequencies, and preparing files for BLAST.

    Parameters
    ----------
    schema_folder : str
        Path to the folder containing schema FASTA files.
    constants : list
        A list of constants used for processing.
    results_output : str
        Path to the directory where results will be saved.

    Returns
    -------
    tuple
        A tuple containing:
        - all_nucleotide_sequences (Dict[str, str]): Dictionary of nucleotide sequences.
        - master_file_path (str): Path to the master FASTA file.
        - trans_paths (Dict[str, str]): Dictionary with the path of the translation file for each loci.
        - to_blast_paths (Dict[str, str]): Dictionary of paths to sequences to be used for BLAST.
        - all_alleles (Dict[str, List[str]]): Dictionary of all alleles with loci IDs as keys.
        - group_reps_ids (Dict[str, List[str]]): Dictionary of group representative IDs.
        - group_alleles_ids (Dict[str, List[str]]): Dictionary of group allele IDs.
        - to_run_against (Dict[str, str]): Dictionary of paths to sequences to be used for BLAST database.
        - new_max_hits (Dict[str, int]): Dictionary with the number of max hits per loci for the BLAST.
        - seqid_file_dict (Dict[str, int]): Dictionary with the paths to the files with the seqid to ignore per loci for the BLAST.
    """
    # Map loci IDs to the paths to the FASTA files in the schema
    loci_files = {file.rsplit('_', 1)[0]: os.path.join(schema_folder, file) for file in os.listdir(schema_folder) if file.endswith('.fasta')}
    # Do the same for the FASTA files containing the representative alleles
    schema_short_dir = os.path.join(schema_folder, 'short')
    representative_files = {file.rsplit('_', 1)[0]: os.path.join(schema_short_dir, file) for file in os.listdir(schema_short_dir) if file.endswith('.fasta')}

    # Path to the master FASTA file
    master_file_path = os.path.join(results_output, 'master_protein.fasta')
    seqid_output = os.path.join(results_output, 'seqid_files')
    ff.create_directory(seqid_output)

    # Create a directory for possible new loci translations
    possible_new_loci_translation_folder = os.path.join(results_output, 'schema_translation_folder')
    ff.create_directory(possible_new_loci_translation_folder)

    # Initialize dictionaries for alleles, translations, and frequencies
    all_alleles: Dict[str, List[str]] = {}
    all_nucleotide_sequences: Dict[str, str] = {} 
    trans_paths: Dict[str, str] = {}
    group_reps_ids: Dict[str, List[str]] = {} 
    group_alleles_ids: Dict[str, List[str]] = {} 
    blast_alleles: Dict[str, List[str]] = {}
    seqid_file_dict: Dict[str, str] = {}

    # Process alleles to run, DNA sequences
    for locus, file in representative_files.items():
        sequences = sf.fetch_fasta_dict(file, False)
        for seqid, sequence in sequences.items():
            group_reps_ids.setdefault(locus, []).append(seqid)
            all_nucleotide_sequences.setdefault(locus, str(sequence))

    # Write master file to run against, DNA sequences
    prf.print_message('Write master file for Blastp.', 'info')

    for loci, loci_path in to_run_against.items():
        negative_seqid_file = os.path.join(seqid_output, f'negative_seqid_{loci_id}.txt')
        blast_alleles.setdefault(loci_id, [])
        seqid_file_dict.setdefault(loci_id, negative_seqid_file)
        fasta_dict = sf.fetch_fasta_dict(loci_path, False)
        for allele_id, sequence in fasta_dict.items():
            blast_alleles[loci_id].append(allele_id)
            group_alleles_ids.setdefault(loci_id, []).append(allele_id)
            all_nucleotide_sequences.setdefault(allele_id, str(sequence))
            protseq = sf.translate_sequence(str(sequence), constants[6])
            # Write to master file
            write_type = 'a' if os.path.exists(master_file_path) else 'w'
            with open(master_file_path, write_type) as master_file:
                master_file.write(f">{allele_id}\n{str(protseq)}\n")
            write_type2 = 'a' if os.path.exists(negative_seqid_file) else 'w'
            with open(negative_seqid_file, write_type2) as seqid_file:
                seqid_file.write(f"{allele_id}\n")

    for loci, loci_path in to_blast_paths.items():
        loci_id = ff.file_basename(loci, False)
        all_alleles.setdefault(loci_id, [])
        fasta_dict = sf.fetch_fasta_dict(loci_path, False)
        for allele_id, sequence in fasta_dict.items():  
            all_alleles[loci_id].append(allele_id)
        # Translate sequences and update translation dictionary
        trans_path_file = os.path.join(possible_new_loci_translation_folder, f"{loci_id}.fasta")
        trans_paths[loci_id] = trans_path_file
        trans_dict, _, _ = sf.translate_seq_deduplicate(fasta_dict,
                                                        trans_path_file,
                                                        constants[5],
                                                        constants[6],
                                                        False)

    return all_nucleotide_sequences, master_file_path, trans_paths, to_blast_paths, all_alleles, group_reps_ids, group_alleles_ids, to_run_against, seqid_file_dict
