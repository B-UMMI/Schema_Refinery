import itertools

def join_list(lst, delimiter):
    """Join all elements in a list into a single string.

    Parameters
    ----------
    lst : list
        List with elements to be joined.
    delimiter : str
        Character used to join list elements.

    Returns
    -------
    joined_list : str
        A single string with all elements in the input
        list joined by the character chosen as link.
    """

    joined_list = delimiter.join(lst)

    return joined_list


def flatten_list(list_to_flatten):
    """Flatten one level of a nested list.

    Parameters
    ----------
    list_to_flatten : list
        Nested list to flatten.

    Returns
    -------
    flattened_list : str
        Input list flattened by one level.
    """

    flattened_list = list(itertools.chain(*list_to_flatten))

    return flattened_list


def isListEmpty(input_list):
    """Check if a nested list is empty."""
    if isinstance(input_list, list):
        return all(map(isListEmpty, input_list)) if isinstance(input_list, list) else False


def divide_list_into_n_chunks(list_to_divide, n):
    """Divides a list into a defined number of sublists.

    Parameters
    ----------
    list_to_divide : list
        List to divide into sublists.
    n : int
        Number of sublists to create.

    Returns
    -------
    sublists : list
        List with the sublists created by dividing
        the input list.
    """
    
    sublists = []
    d, r = divmod(len(list_to_divide), n)
    for i in range(n):
        si = (d+1)*(i if i < r else r) + d*(0 if i < r else i - r)
        sublists.append(list_to_divide[si:si+(d+1 if i < r else d)])

    # exclude lists that are empty due to small number of elements
    sublists = [line for line in sublists if len(line) > 0]

    return sublists

def get_max_min_values(input_list):
    """
    From an input list return the max and min int/float present

    Parameters
    ----------
    input_list: int
        List containing values
    
    Returns
    -------
    Returns largest value and smallest value
    """

    return max(input_list), min(input_list)