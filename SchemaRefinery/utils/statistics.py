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

def two_sets_within_min_max_threshold(set1: List[float], set2: List[float], threshold: float) -> Tuple[bool, bool]:
    """
    Check if the minimum and maximum values of two sets are within a specified threshold of each other.

    Parameters
    ----------
    set1 : list of float
        The first set of min and max values.
    set2 : list of float
        The second set of min and max values.
    threshold : float
        The threshold to compare against (e.g., 0.2 for 20%).

    Returns
    -------
    tuple of bool
        A tuple containing two boolean values:
        - The first boolean indicates if the minimum values of the two sets are within the specified threshold.
        - The second boolean indicates if the maximum values of the two sets are within the specified threshold.

    Examples
    --------
    >>> two_sets_within_min_max_threshold([10, 20, 30], [12, 22, 32], 0.2)
    (True, True)
    >>> two_sets_within_min_max_threshold([10, 20, 30], [15, 25, 35], 0.2)
    (False, True)
    """

    min_within_x_threshold = abs(set1[0] - set2[0]) <= threshold * max(set1[0], set2[0])
    max_within_x_threshold = abs(set1[1] - set2[1]) <= threshold * max(set1[1], set2[1])


    return min_within_x_threshold, max_within_x_threshold

