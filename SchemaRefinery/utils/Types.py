from typing import Dict, List, TypedDict, Tuple, Set, Any, Union
from collections import OrderedDict

MergedAllClasses = Dict[str, Union[List[str], Dict[str, List[str]], Set[str]]]

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
    classification: str 

BlastDict = Dict[str, Dict[str, List[BlastResult]]]

RepresentativeBlastResultsCoords = Dict[str, Dict[str, Dict[str, List[Tuple[int, int]]]]]

BSRValues = Dict[str, Dict[str, float]]

class ProcessedResult(TypedDict):
    class_: str
    ids: List[str]
    empty_list: List[str]
    pair: Tuple[str, str]
    combined_list: List[str]

ProcessedResults = Dict[str, Tuple[ProcessedResult]]

CountResultsByClass = Dict[str, Dict[str, OrderedDict[str, int]]]

CountResultsByClassWithInverse = Dict[str, Dict[str, Dict[str, List[Dict[str, int]]]]]

RepsAndAllelesIds = Dict[str, Tuple[Set[str], Set[str]]]

AllRelationships = Dict[str, List[List[str]]]

RelatedClusters = Dict[str, Dict[int, List[str]]]

Recomendations = Dict[str, Dict[int, Dict[str, Set[str]]]]

class Metadata(TypedDict):
    total_count: int
    reports: List[Dict[str, Any]]