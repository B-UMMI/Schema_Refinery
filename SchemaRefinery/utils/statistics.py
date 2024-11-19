from typing import List, Tuple

def modes_within_value(mode1: float, mode2: float, value: float) -> bool:
    """
    Check if the absolute difference between two modes is within a specified value.

    Parameters
    ----------
    mode1 : float
        The first mode value.
    mode2 : float
        The second mode value.
    value : float
        The value to compare against.

    Returns
    -------
    bool
        True if the absolute difference between mode1 and mode2 is less than or equal to the specified value times the maximum of mode1 and mode2, False otherwise.

    Examples
    --------
    >>> modes_within_value(10, 12, 0.2)
    True
    >>> modes_within_value(10, 15, 0.2)
    False
    """
    return abs(mode1 - mode2) <= value * max(mode1, mode2)

def if_loci_intersect(set1: List[float], set2: List[float]) -> bool:
    """
    Check if two loci (intervals) intersect.

    Parameters
    ----------
    set1 : list of float
        The first set of values representing a locus (interval).
    set2 : list of float
        The second set of values representing a locus (interval).

    Returns
    -------
    bool
        True if the loci intersect, False otherwise.

    Examples
    --------
    >>> if_loci_intersect([1.0, 5.0], [4.0, 6.0])
    True
    >>> if_loci_intersect([1.0, 3.0], [4.0, 6.0])
    False
    """
    # Check if the intervals intersect
    return set1[0] <= set2[1] and set1[1] >= set2[0] or set2[0] <= set1[1] and set2[1] >= set1[0]

def calculate_loci_distance(set1: List[float], set2: List[float], threshold: float) -> bool:
    """
    Calculate if the distance between two loci (intervals) is within a specified threshold.

    Parameters
    ----------
    set1 : list of float
        The first set of values representing a locus (interval).
    set2 : list of float
        The second set of values representing a locus (interval).
    threshold : float
        The threshold to compare against.

    Returns
    -------
    bool
        True if the distance between the loci is within the specified threshold, False otherwise.

    Examples
    --------
    >>> calculate_loci_distance([1.0, 5.0, 2.0], [4.0, 6.0, 3.0], 0.2)
    True
    >>> calculate_loci_distance([1.0, 3.0, 2.0], [4.0, 6.0, 3.0], 0.2)
    False
    """
    # Calculate the minimum distance between the intervals and check if it is within the threshold
    return min(abs(set1[0] - set2[1]), abs(set2[0] - set1[1])) <= max(set1[2], set2[2]) * threshold