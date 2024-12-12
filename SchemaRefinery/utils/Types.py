from typing import Dict, List, TypedDict, Tuple

class ClustersToKeep(TypedDict):
    '1a': Dict[int, List[str]]
    '1b': List[str]
    '2a': List[str]
    '3a': List[str]
    '2b': List[str]
    '1c': List[str]
    '3b': List[str]
    '4a': List[str]
    '4b': List[str]
    '4c': List[str]
    '5': List[str]
    'Retained_not_matched_by_blastn': List[str]

class BlastResult(TypedDict):
    query: str
    subject: str
    query_length: int
    subject_length: int
    query_start: int
    query_end: int
    subject_start: int
    subject_end: int
    length: int
    score: int
    gaps: int
    pident: float
    bsr: float
    kmers_sim: int
    kmers_cov: int
    frequency_in_genomes_query_cds: int
    frequency_in_genomes_subject_cds: int
    global_palign_all_min: float
    global_palign_all_max: float
    global_palign_pident_min: float
    global_palign_pident_max: float
    local_palign_min: float
    class_: str 

class BlastDict(TypedDict):
    query: Dict[str, Dict[int, BlastResult]]

class RepresentativeBlastResultsCoords(TypedDict):
    query: Dict[str, Dict[str, List[Tuple[int, int]]]]
