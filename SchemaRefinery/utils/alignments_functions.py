from copy import deepcopy
from typing import List, Tuple, Dict, Any, Set, Union
try:
    from utils import (iterable_functions as itf,
                       constants as ct,
                       Types as tp,
                       print_functions as pf)
except ModuleNotFoundError:
    from SchemaRefinery.utils import (iterable_functions as itf,
                                      constants as ct,
                                      Types as tp,
                                      print_functions as pf)

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

            if second["start"] - first["stop"] <= ct.MAX_GAP_UNITS:
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
                found_new_interval = True
                start_stop_list_for_processing = start_stop_list_for_processing[:i] + [new_interval] + start_stop_list_for_processing[i + 2:]
                new_index = i
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
                has_alignment = True
                break
        
        if not has_alignment:
            new_original.append(i)

    return new_original


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


def get_alignments_dict_from_blast_results(
    blast_results_file: str,
    pident_threshold: float,
    get_coords: bool,
    get_self_score: bool,
    skip_reverse_alignments: bool,
    if_alleles: bool,
    multiple_reps: bool
) -> Tuple[tp.BlastDict, Union[int, Dict[str, int]], tp.RepresentativeBlastResultsCoords, tp.RepresentativeBlastResultsCoords]:
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

    alignments_dict: tp.BlastDict = {}
    alignment_coords_pident: tp.RepresentativeBlastResultsCoords = {}
    alignment_coords_all: tp.RepresentativeBlastResultsCoords = {}
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


def get_alignments_dict_from_blast_results_simplified(
    blast_results_file: str,
    pident_threshold: float,
    get_coords: bool,
    skip_reverse_alignments: bool
) -> Tuple[tp.BlastDict, tp.RepresentativeBlastResultsCoords, tp.RepresentativeBlastResultsCoords]:
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
    skip_reverse_alignments : bool
        Whether to skip inverse alignments.

    Returns
    -------
    Tuple[tp.BlastDict, tp.RepresentativeBlastResultsCoords, tp.RepresentativeBlastResultsCoords]
        - alignments_dict: Dictionary containing the results of the BLAST.
        - alignment_coords_all: Contains the coordinates for the query/subject pair, the coordinates are in reference to the query.
        - alignment_coords_pident: Contains the coordinates for the query/subject pair, the coordinates are in reference to the query and filtered by pident threshold.
    """

    alignments_dict: tp.BlastDict = {}
    alignment_coords_pident: tp.RepresentativeBlastResultsCoords = {}
    alignment_coords_all: tp.RepresentativeBlastResultsCoords = {}

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

    return alignments_dict, alignment_coords_all, alignment_coords_pident


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