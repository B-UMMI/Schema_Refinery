from copy import deepcopy
from typing import List, Tuple, Dict, Any, Union
try:
    from utils import (iterable_functions as itf,
                       Types as tp,)
except ModuleNotFoundError:
    from SchemaRefinery.utils import (iterable_functions as itf,
                                      Types as tp,)

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