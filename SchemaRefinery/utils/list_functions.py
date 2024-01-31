import itertools
from collections import defaultdict

def join_list(lst, delimiter):
    """
    Join all elements in a list into a single string.

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
    """
    Flatten one level of a nested list.

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
    """
    Check if a nested list is empty.
    """
    
    if isinstance(input_list, list):
        return all(map(isListEmpty, input_list)) if isinstance(input_list, list) else False


def divide_list_into_n_chunks(list_to_divide, n):
    """
    Divides a list into a defined number of sublists.

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

def decompress_number(text, index):
    """
    Decode a single number from a string created with polyline encoding.

    Parameters
    ----------
    text : str
        String representing a compressed list of numbers.
    index : int
        Index of the first character to start decoding a number.

    Returns
    -------
    Index to start decoding the next number and the number decoded
    in the current function call.
    """
    number = 0
    bitwise_shift = 0

    while True:
        # subtract 63 and remove bit from OR with 0x20 if 5-bit chunk
        # is not the last to decode a number that was in the original list
        n = (ord(text[index]) - 63)
        index += 1
        # only continue if there is a continuation bit (0b100000)
        # e.g.: 0b11111 only has 5 bits, meaning it is the last chunk
        # that is necessary to decode the number
        if n >= 0x20:
            # subtract 0x20 (0b100000) to get the 5-bit chunk
            n -= 0x20
            # contruct the binary number with biwise shift to add each
            # 5-bit chunk to original position
            number = number | (n << bitwise_shift)
            # increment bitwise shift value for next 5-bit chunk
            bitwise_shift += 5
        else:
            break

    # add the last chunk, without continuation bit, to the leftmost position
    number = number | (n << bitwise_shift)

    # invert bits to get negative number if sign bit is 1
    # remove sign bit and keep only bits for decoded value
    return index, (~number >> 1) if (number & 1) != 0 else (number >> 1)

def polyline_decoding(text):
    """
    Decode a list of integers compressed with polyline encoding.

    Parameters
    ----------
    text : str
        String representing a compressed list of numbers.

    Returns
    -------
    number_list : list
        List with the decoded numbers.
    """
    
    number_list = []
    index = last_num = 0
    # decode precision value
    precision = ord(text[index]) - 63
    index += 1
    while index < len(text):
        # decode a number and get index to start decoding next number
        index, difference = decompress_number(text, index)
        # add decoded difference to get next number in original list
        last_num += difference
        number_list.append(last_num)

    number_list = [round(item * (10 ** (-precision)), precision) for item in number_list]

    return number_list

def merge_dicts(dicts_list):
    """
    Function that merges a dict of dicts into one single dict with only the first
    level.
    
    Parameters
    ----------
    dicts_list : dict
        Dict of dicts.
        
    Returns
    -------
    out_dict : dict
        Dict with merged entries from input dict at the first level of the dict.
    """
    
    out_dict = defaultdict(dict)
    for dict_ in dicts_list: # you can list as many input dicts as you want here
        for key, value in dicts_list[dict_].items():
            out_dict[key].update(value)
            
    return out_dict