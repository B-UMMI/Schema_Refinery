from copy import deepcopy
from typing import List, Tuple, Dict, Any, Set, Union
try:
    from RefineSchema.constants import MAX_GAP_UNITS
    from utils import iterable_functions as itf
except ModuleNotFoundError:
    from SchemaRefinery.RefineSchema.constants import MAX_GAP_UNITS
    from SchemaRefinery.utils import iterable_functions as itf

def join_intervals(alignments: List[List[Any]]) -> Tuple[List[str], List[Dict[str, Any]]]:
    """
    Joins alignments that intersect with each other and merges them into a single alignment,
    taking into account custom scoring that needs to be calculated for the alleles.

    Parameters
    ----------
    alignments : List[List[Any]]
        List containing the data obtained from BLAST. Each alignment is expected to be a list with the following elements:
        [start, stop, ..., pident, gaps, length, query, subject, internal_alignment]

    Returns
    -------
    Tuple[List[str], List[Dict[str, Any]]]
        A tuple containing:
        - final_start_stop_list: List of strings representing the final start-stop intervals with joined intervals.
        - start_stop_list_for_processing: List of dictionaries containing detailed information about each interval.
    """
    
    start_stop_list_for_processing: List[Dict[str, Any]] = [
        {
            "start": alignment[0], 
            "stop": alignment[1], 
            "joined_intervals": set(), 
            "pident": alignment[4] / 100, 
            "gaps": alignment[5], 
            "length": alignment[6],
            "product_length_pident_list": [alignment[4] / 100 * alignment[6]], 
            "length_list": [alignment[6]],
            "query": alignment[7],
            "subject": alignment[8],
            "internal_alignments": [alignment[9]]
        } 
        for alignment in alignments
    ]
    
    found_new_interval: bool = True
    new_index: int = 0
    while found_new_interval:
        found_new_interval = False
        i: int
        for i in range(new_index, len(start_stop_list_for_processing) - 1):
            first = start_stop_list_for_processing[i]
            second = start_stop_list_for_processing[i + 1]

            if second["start"] - first["stop"] <= MAX_GAP_UNITS:
                new_last: int = max(first["stop"], second["stop"])
                new_first: int = min(first["start"], second["start"])

                new_joined_intervals: Set[Tuple[int, int]] = first["joined_intervals"].union(
                    second["joined_intervals"]
                ).union({(first['start'], first['stop']), (second['start'], second['stop'])})
                
                new_interval: Dict[str, Any] = {
                    "start": new_first, 
                    "stop": new_last, 
                    "joined_intervals": new_joined_intervals, 
                    "pident": f"{first['pident']};{second['pident']}", 
                    "gaps": f"{first['gaps']};{second['gaps']}",
                    "product_length_pident_list": first["product_length_pident_list"] + second["product_length_pident_list"], 
                    "length_list": first["length_list"] + second["length_list"],
                    "query": first["query"],
                    "subject": first["subject"],
                    "internal_alignments": first["internal_alignments"] + second["internal_alignments"]
                }
                found_new_interval: bool = True
                start_stop_list_for_processing: List[Dict[str, Any]]
                new_interval: Dict[str, Any]
                start_stop_list_for_processing = start_stop_list_for_processing[:i] + [new_interval] + start_stop_list_for_processing[i + 2:]
                new_index: int = i
                break

    for interval in start_stop_list_for_processing:
        interval["custom_scoring"] = sum(interval["product_length_pident_list"]) / sum(interval["length_list"])
        del interval["product_length_pident_list"]
        interval['joined_intervals'] = sorted(list(interval['joined_intervals']), key=lambda x: x[0])
        interval['joined_intervals'] = [f"{i[0]}-{i[1]}" for i in interval['joined_intervals']]
        interval['joined_intervals'] = f"({';'.join(interval['joined_intervals'])})" if interval['joined_intervals'] else ""

    final_start_stop_list: List[str] = [
        f"{interval['start']}-{interval['stop']}{interval['joined_intervals']}" 
        for interval in start_stop_list_for_processing
    ]
    
    return final_start_stop_list, start_stop_list_for_processing

def filter_out_equal_alignments(original: List[Dict[str, Any]], inverted: List[Dict[str, Any]]) -> List[Dict[str, Any]]:
    """
    Filters out alignments from the inverted list that are equal to any alignment in the original list.
    An alignment is considered equal if it has the same start and stop positions but with the query and subject inverted.

    Parameters
    ----------
    original : List[Dict[str, Any]]
        List containing alignments related to one key.
    inverted : List[Dict[str, Any]]
        List containing alignments related to the inverse of the previous key.

    Returns
    -------
    List[Dict[str, Any]]
        A new list containing the original alignments with the non-equal inverted alignments appended.
    """
    
    new_original: List[Dict[str, Any]] = deepcopy(original)

    for i in inverted:
        has_alignment: bool = False
        for o in original:
            if i["query_start"] == o["subject_start"] and i["query_end"] == o["subject_end"]:
                has_alignment: bool = True
                break
        
        if not has_alignment:
            new_original.append(i)

    return new_original

def process_alignments_for_graphs(alignments_dict: Dict[str, List[Dict[str, Any]]]) -> Dict[str, List[Dict[str, Any]]]:
    """
    Processes alignments to join alignments from an inverse of another alignment, if it exists, so they appear in the same graph.

    Parameters
    ----------
    alignments_dict : Dict[str, List[Dict[str, Any]]]
        Dictionary containing alignments. The keys are strings representing loci pairs separated by a semicolon,
        and the values are lists of dictionaries containing alignment details.

    Returns
    -------
    Dict[str, List[Dict[str, Any]]]
        Processed dictionary with alignments joined where inverses exist.
    """

    keys_set: set = set()
    processed_alignments_dict: Dict[str, List[Dict[str, Any]]] = {}

    print("Processing alignments for graphs...")

    for key in alignments_dict.keys():
        locus_1, locus_2 = key.split(";")
        inverted_key = f"{locus_2};{locus_1}"

        if inverted_key in alignments_dict.keys():
            if key in keys_set:
                continue
            else:
                original_alignments: dict = alignments_dict[key]
                inverted_alignments: dict = alignments_dict[inverted_key]

                filtered_alignments: dict = filter_out_equal_alignments(original_alignments, inverted_alignments)

                processed_alignments_dict[key] = filtered_alignments

        else:
            processed_alignments_dict[key] = alignments_dict[key]

        keys_set.add(inverted_key)
    
    return processed_alignments_dict

def get_alignments_dict(blast_results_file: str) -> Dict[str, List[Dict[str, Any]]]:
    """
    Organizes alignments with the same key "Locus_A:Locus_B" into the same dictionary.
    Builds a dictionary where the key is "Locus_A:Locus_B" and the value is a list of all alignments for that key correspondence.

    Parameters
    ----------
    blast_results_file : str
        Path to the BLAST result file.

    Returns
    -------
    Dict[str, List[Dict[str, Any]]]
        Dictionary containing the necessary information for graph building and representatives vs alleles BLAST.
    """

    alignments_dict: Dict[str, List[Dict[str, Any]]] = {}
    with open(blast_results_file, "r") as f:
        lines = f.readlines()
        for line in lines:
            cols: List[str] = line.strip().split("\t")
            query: str = cols[0]
            subject: str = cols[1]
            query_length: int = int(cols[2])
            subject_length: int = int(cols[3])
            query_start: int = int(cols[4])
            query_end: int = int(cols[5])
            subject_start: int = int(cols[6])
            subject_end: int = int(cols[7])
            length: int = int(cols[8])
            score: int = int(cols[9])
            gaps: int = int(cols[10])
            pident: float = float(cols[11])

            key: str = f"{query};{subject}"
            value: dict = {
                "query": query,
                "subject": subject,
                "query_length": int(query_length),
                "subject_length": int(subject_length),
                "query_start": int(query_start),
                "query_end": int(query_end),
                "subject_start": int(subject_start),
                "subject_end": int(subject_end),
                "length": int(length),
                "score": int(score),
                "gaps": int(gaps),
                "pident": float(pident)
            }

            if key not in alignments_dict:
                alignments_dict[key] = [value]
            else:
                alignments_dict[key].append(value)

    return alignments_dict

def process_blast_results(blast_results_file: str, constants_threshold: List[float]) -> Tuple[str, Dict[str, List[Dict[str, Any]]]]:
    """
    Main function to process the received BLAST results. Filters the results, organizes the alignments, and returns 
    a string with the information of the selected alignments for the report and a dictionary with all the selected 
    alignments to build the graphs.

    Parameters
    ----------
    blast_results_file : str
        Path to the BLAST result file.
    constants_threshold : List[float]
        List that contains two constants: alignment_ratio_threshold and pident_threshold.

    Returns
    -------
    Tuple[str, Dict[str, List[Dict[str, Any]]]]
        A tuple containing:
        - alignment_strings: str
            Used to write inside the alleles_report file.
        - filtered_alignments_dict: Dict[str, List[Dict[str, Any]]]
            Used downstream to generate graphs.
    """

    alignments_dict: Dict[str, List[Dict[str, Any]]] = get_alignments_dict(blast_results_file)

    alignment_ratio_threshold: float
    pident_threshold: float
    alignment_ratio_threshold, pident_threshold = constants_threshold

    # Filter alignments by pident
    alignments_dict: dict  = {
        key: [alignment for alignment in alignments if alignment["pident"] >= pident_threshold]
        for key, alignments in alignments_dict.items()
    }
    # Remove dictionary entries with zero alignments after filtering by pident
    alignments_dict: dict = {
        key: alignments for key, alignments in alignments_dict.items() if len(alignments) != 0
    }

    filtered_alignments_dict: Dict[str, List[Dict[str, Any]]] = deepcopy(alignments_dict)
    alignment_strings: List[str] = []

    for key, alignments in alignments_dict.items():
        if len(alignments) > 0:
            query: str = alignments[0]['query']
            subject: str = alignments[0]['subject']
            
            query_before_underscore: str = query.split("_")[0]
            subject_before_underscore: str = subject.split("_")[0]

            if query_before_underscore == subject_before_underscore:
                del filtered_alignments_dict[key]
                continue
            else:
                query_length = alignments[0]["query_length"]
                subject_length = alignments[0]["subject_length"]

                alignments.sort(key=lambda x: x["query_start"])
                query_start_stops_list = [
                    [entry["query_start"], entry["query_end"], entry["query_length"], entry["subject_length"], entry["pident"], entry["gaps"], entry["length"], entry["query"], entry["subject"], entry]
                    for entry in alignments
                ]
                alignments.sort(key=lambda x: x["subject_start"])
                subject_start_stops_list = [
                    [entry["subject_start"], entry["subject_end"], entry["query_length"], entry["subject_length"], entry["pident"], entry["gaps"], entry["length"], entry["query"], entry["subject"], entry]
                    for entry in alignments
                ]

                final_query_start_stop_list, alignment_query = join_intervals(query_start_stops_list)
                final_subject_start_stop_list, alignment_subject = join_intervals(subject_start_stops_list)
                
                final_gaps = ';'.join([str(alignment["gaps"]) for alignment in alignment_query])
                final_pident = ';'.join([str(alignment["pident"]) for alignment in alignment_query])

                alignment_query.sort(key=lambda x: (x["stop"] - x["start"]), reverse=True)
                alignment_subject.sort(key=lambda x: (x["stop"] - x["start"]), reverse=True)
                bigger_query_alignment = alignment_query[0]["stop"] - alignment_query[0]["start"]
                bigger_subject_alignment = alignment_subject[0]["stop"] - alignment_subject[0]["start"]

                query_ratio: float = bigger_query_alignment / query_length
                subject_ratio: float = bigger_subject_alignment / subject_length

                if query_ratio >= alignment_ratio_threshold or subject_ratio >= alignment_ratio_threshold:
                    query_start_stops: str = ';'.join(final_query_start_stop_list)
                    subject_start_stops: str = ';'.join(final_subject_start_stop_list)
                    alignment_string: str = f"{query}\t{subject}\t{query_start_stops}\t{subject_start_stops}\t{query_ratio}\t{subject_ratio}\t{query_length}\t{subject_length}\t{final_gaps}\t{final_pident}\n"
                    alignment_strings.append(alignment_string)
                else:
                    del filtered_alignments_dict[key]
        
    return '\n'.join(alignment_strings), filtered_alignments_dict

def get_alignments_dict_from_blast_results(
    blast_results_file: str,
    pident_threshold: float,
    get_coords: bool,
    get_self_score: bool,
    skip_reverse_alignments: bool,
    if_alleles: bool,
    multiple_reps: bool
) -> Tuple[Dict[str, Dict[str, Dict[int, Dict[str, Any]]]], Union[int, Dict[str, int]], Dict[str, Dict[str, Dict[str, List[List[int]]]]], Dict[str, Dict[str, Dict[str, List[List[int]]]]]]:
    """
    Reads BLAST results file and extracts the necessary items. Based on input, also fetches the coordinates 
    based on query sequences and self-score contained inside the BLAST results file.

    Parameters
    ----------
    blast_results_file : str
        Path to the BLAST result file.
    pident_threshold : float
        Pident threshold to exclude BLAST results.
    get_coords : bool
        Whether to fetch coordinates for the BLAST match.
    get_self_score : bool
        Whether to get self-score from BLAST results (Note: for self-score to be fetched, the query must also be in the subjects database).
    skip_reverse_alignments : bool
        Whether to skip inverse alignments.
    if_alleles : bool
        If True, the function will process alleles to fetch self-score.
    multiple_reps : bool
        If True, the function will process multiple representatives to fetch self-score; otherwise, it will fetch the highest self-score.

    Returns
    -------
    Tuple[Dict[str, Dict[str, Dict[int, Dict[str, Any]]]], Union[int, Dict[str, int]], Dict[str, Dict[str, Dict[str, List[List[int]]]]], Dict[str, Dict[str, Dict[str, List[List[int]]]]]]
        - alignments_dict: Dictionary containing the results of the BLAST.
        - self_score: Value of self-score inside the BLAST results file.
        - alignment_coords_all: Contains the coordinates for the query/subject pair, the coordinates are in reference to the query.
        - alignment_coords_pident: Contains the coordinates for the query/subject pair, the coordinates are in reference to the query and filtered by pident threshold.

    Notes
    -----
    The alignments dicts (alignments_dict, alignment_coords_all, alignment_coords_pident) will have the following structure:
    alignments_dict = {query: {subject: {1: {alignment1}, 2: {alignment2}, ...}, ...}, ...}
    where query will be the allele identifier e.g., loci1_1, so for the same loci there may be various query entries.
    The self_score will be the highest self-score found in the BLAST results file or a dictionary containing the self-score for each representative present.
    """

    alignments_dict: Dict[str, Dict[str, Dict[int, Dict[str, Any]]]] = {}
    alignment_coords_pident: Dict[str, Dict[str, Dict[str, List[List[int]]]]] = {}
    alignment_coords_all: Dict[str, Dict[str, Dict[str, List[List[int]]]]] = {}
    pattern: str = '_(\d+)'
    self_scores: Union[int, Dict[str, int]] = 0 if not multiple_reps else {}

    with open(blast_results_file, "r") as f:
        lines: List[str] = f.readlines()
        i: int = 1
        for line in lines:
            # Extract the columns into the variables
            cols: List[str] = line.strip().split("\t")
            query: str = cols[0]
            subject: str = cols[1]
            query_length: int = int(cols[2])
            subject_length: int = int(cols[3])
            query_start: int = int(cols[4])
            query_end: int = int(cols[5])
            subject_start: int = int(cols[6])
            subject_end: int = int(cols[7])
            length: int = int(cols[8])
            score: int = int(cols[9])
            gaps: int = int(cols[10])
            pident: float = float(cols[11])

            # Save the dict
            value: Dict[str, Any] = {
                "query": query,
                "subject": subject,
                "query_length": query_length,
                "subject_length": subject_length,
                "query_start": query_start,
                "query_end": query_end,
                "subject_start": subject_start,
                "subject_end": subject_end,
                "length": length,
                "score": score,
                "gaps": gaps,
                "pident": pident
            }
            
            if if_alleles:
                if itf.remove_by_regex(query, pattern) == itf.remove_by_regex(subject, pattern):
                    if not multiple_reps:
                        self_score = 0
                    else:
                        self_scores.setdefault(query, 0)
                        self_score = self_scores[query]
                    # Largest self-score is chosen
                    if pident == 100 and get_self_score and score > self_score:
                        if multiple_reps:
                            self_scores[query] = score
                        else:
                            self_scores = score
                    continue
            # Skip if entry matched itself and get self-score if needed
            elif query == subject:
                if not multiple_reps:
                    self_score = 0
                else:
                    self_scores.setdefault(query, 0)
                    self_score = self_scores[query]
                # Largest self-score is chosen
                if pident == 100 and get_self_score and score > self_score:
                    if multiple_reps:
                        self_scores[query] = score
                    else:
                        self_scores = score
                continue
            
            if skip_reverse_alignments:
                if query_start > query_end or subject_start > subject_end:
                    continue

            if query not in alignments_dict:
                alignments_dict[query] = {}
                if get_coords:
                    alignment_coords_all[query] = {}
                    alignment_coords_pident[query] = {}
            if subject not in alignments_dict[query]:
                # Create and save the first entry of BLAST
                alignments_dict[query][subject] = {i: value}
                if get_coords:
                    alignment_coords_all[query][subject] = {
                        'query': [[query_start, query_end]], 
                        'subject': [[subject_start, subject_end]]
                    }
                    # Align by pident
                    if pident >= pident_threshold:
                        alignment_coords_pident[query][subject] = {
                            'query': [[query_start, query_end]],
                            'subject': [[subject_start, subject_end]]
                        }
                    else:
                        alignment_coords_pident[query][subject] = {
                            'query': [],
                            'subject': []
                        }
            else:
                # Save the other entries based on total number of entries present to get the ID
                k: int = max(alignments_dict[query][subject].keys()) + 1
                alignments_dict[query][subject][k] = value
                if get_coords:
                    alignment_coords_all[query][subject]['query'].append([query_start, query_end])
                    alignment_coords_all[query][subject]['subject'].append([subject_start, subject_end])
                    # Align by pident
                    if pident >= pident_threshold:
                        alignment_coords_pident[query][subject]['query'].append([query_start, query_end])
                        alignment_coords_pident[query][subject]['subject'].append([subject_start, subject_end])
            
    return alignments_dict, self_scores, alignment_coords_all, alignment_coords_pident

def remove_inverse_alignments(
    alignments_dict: Dict[str, Dict[str, Any]], 
    all_representatives_alignments_dict: Dict[str, Dict[str, Any]]
) -> List[List[str]]:
    """
    Filters out inverse alignments from the alignments dictionary to avoid repeated results in the alleles report.

    Parameters
    ----------
    alignments_dict : Dict[str, Dict[str, Any]]
        Dictionary containing alignments of the representative loci.
    all_representatives_alignments_dict : Dict[str, Dict[str, Any]]
        Dictionary containing all representative loci alignments without their inverses.

    Returns
    -------
    List[List[str]]
        List of lists, where each sublist contains the loci alignment pair IDs.
    """

    filtered_alignments_dict: Dict[str, Dict[str, Any]] = deepcopy(alignments_dict)
    for key in alignments_dict.keys():
        query: str
        subject: str
        query, subject = key.split(";")
        inverse_key: str = f"{subject};{query}"
        if inverse_key in all_representatives_alignments_dict:
            del filtered_alignments_dict[key]

    all_representatives_alignments_dict.update(filtered_alignments_dict)

    alignments_pair_list: List[List[str]] = [key.split(";") for key in filtered_alignments_dict.keys()]

    return alignments_pair_list

def merge_intervals(intervals: List[List[int]]) -> List[List[int]]:
    """
    Merges intersecting intervals.

    Parameters
    ----------
    intervals : List[List[int]]
        List of intervals, where each interval is represented as a list of two integers [start, end].

    Returns
    -------
    List[List[int]]
        List of merged intervals, where overlapping intervals are combined into a single interval.
    """

    if not intervals:
        return []

    merged: List[List[int]] = [deepcopy(intervals[0])]
    for current in intervals[1:]:
        previous: List[int] = merged[-1]
        # current and previous intervals intersect
        if current[0] <= previous[1]:
            # determine top position
            previous[1] = max(previous[1], current[1])
        # current and previous intervals do not intersect
        else:
            merged.append(deepcopy(current))

    return merged