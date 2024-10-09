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