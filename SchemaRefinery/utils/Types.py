from typing import Dict, List, TypedDict

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