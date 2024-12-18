from typing import Dict, List, TypedDict, Tuple, Set, Any
from collections import OrderedDict

class ClustersToKeep(TypedDict, total=False):
    a1: Dict[int, List[str]]
    b1: List[str]
    a2: List[str]
    a3: List[str]
    b2: List[str]
    c1: List[str]
    b3: List[str]
    a4: List[str]
    b4: List[str]
    c4: List[str]
    five: List[str]
    retained_not_matched_by_blastn: List[str]

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

class BSRValues(TypedDict):
    query_id: Dict[str, float]

class ProcessedResult(TypedDict):
    class_: str
    ids: List[str]
    empty_list: List[str]
    pair: Tuple[str, str]
    combined_list: List[str]

class ProcessedResults(TypedDict):
    query_subject: Tuple[ProcessedResult]

class CountResultsByClass(TypedDict):
    query_subject: Dict[str, OrderedDict[str, int]]

class ClassCount(TypedDict):
    direct_class: int
    inverse_class: int

class CountResultsByClassWithInverse(TypedDict):
    query_subject: Dict[str, Dict[str, List[ClassCount]]]

class RepsAndAllelesIds(TypedDict):
    query_subject: Tuple[Set[str], Set[str]]

class AllRelationships(TypedDict):
    class_: List[List[str]]

class RelatedClusters(TypedDict):
    cluster_id: Dict[int, List[str]]

class Recomendations(TypedDict):
    cluster_id: Dict[int, Dict[str, Set[str]]]

class Metadata(TypedDict):
    total_count: int
    reports: List[Dict[str, Any]]