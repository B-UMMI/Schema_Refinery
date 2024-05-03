import itertools
import pickle

def join_list(lst, delimiter):
    """
    Join all elements in a list into a single string using a delimiter.

    Parameters
    ----------
    lst : list
        List with elements to be joined.
    delimiter : str
        Character or string used to join list elements.

    Returns
    -------
    joined_list : str
        A single string with all elements in the input list joined by the specified delimiter.
    """

    joined_list = delimiter.join(map(str, lst))

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
    
    Parameters
    ----------
    input_list : list
        Input list to verify
    Returns
    -------
    return : bool
        If list is empty or not.
    """
    
    if isinstance(input_list, list):
        return all(map(isListEmpty, input_list)) if isinstance(input_list, list) else False


def divide_list_into_n_chunks(list_to_divide, n):
    """
    Divides a list into a specified number of sublists.

    Parameters
    ----------
    list_to_divide : list
        The list to divide into sublists.
    n : int
        The number of sublists to create.

    Returns
    -------
    sublists : list
        A list of sublists created by dividing the input list.
    """
    
    sublists = []
    list_length = len(list_to_divide)
    sublist_size, remainder = divmod(list_length, n)
    
    # Determine sublist sizes and create sublists
    start_idx = 0
    for i in range(n):
        sublist_length = sublist_size + (1 if i < remainder else 0)
        end_idx = start_idx + sublist_length
        sublists.append(list_to_divide[start_idx:end_idx])
        start_idx = end_idx

    return [sublist for sublist in sublists if sublist]  # Filter out empty sublists

def get_max_min_values(input_list):
    """
    From an input list, return the maximum and minimum integer or float values.

    Parameters
    ----------
    input_list : list[int] or list[float]
        List containing integer or float values.

    Returns
    -------
    return : tuple
        A tuple containing the largest and smallest values in the input list.
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
    return : int
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
            # contruct the binary number with bitwise shift to add each
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

def decode_CDS_sequences_ids(path_to_file):
    """
    Read a dictionary contained in a pickle file and decode its values based on polyline.

    Parameters
    ----------
    path_to_file : str
        Path to the pickle file.

    Returns
    -------
    decoded_dict : dict
        Contains hashed CDS as keys, and number id of the genome where that CDS is present.
    """
    # Load pickle file
    with open(path_to_file, "rb") as infile:
        hash_table = pickle.load(infile)

    # Decode the values using polyline
    decoded_dict = {key: polyline_decoding(value) for key, value in hash_table.items()}

    return decoded_dict

def get_unique_sublists(list_of_lists):
    """
    Identify unique sublists within a list of lists.

    Parameters
    ----------
    list_of_lists : list
        The list containing various sublists.

    Returns
    -------
    unique_sublists : list
        List containing unique sublists.
    """
    seen = set()
    unique_sublists = []
    for sublist in list_of_lists:
        sublist_tuple = tuple(sublist)
        if sublist_tuple not in seen:
            seen.add(sublist_tuple)
            unique_sublists.append(sublist)
    return unique_sublists

def all_match_lists(list1, list2):
    """
    Check if all elements of list1 are contained in list2.

    Parameters
    ----------
    list1 : list
        The query list to check for full containment.
    list2 : list
        The subject list to compare against the query list.

    Returns
    -------
    return : bool
        True if all elements of list1 are found in list2, False otherwise.
    """
    return all(elem in list2 for elem in list1)

def any_match_lists(list1, list2):
    """
    Check if any element of list1 is contained in list2.

    Parameters
    ----------
    list1 : list
        The query list to check for partial containment.
    list2 : list
        The subject list to compare against the query list.

    Returns
    -------
    return : bool
        True if any element of list1 is found in list2, False otherwise.
    """
    return any(elem in list2 for elem in list1)

def contains_sublist(main_list, list_of_lists):
    """
    Check if elements of the main list are fully contained in any sublist of the list of lists.

    Parameters
    ----------
    main_list : list
        The query list to check for full containment.
    list_of_lists : list
        The subject list containing sublists to compare against the query list.

    Returns
    -------
    return : bool
        True if all elements of the main list are fully contained in any sublist of the list of lists, False otherwise.
    """
    for sublist in list_of_lists:
        if all_match_lists(main_list, sublist):
            return True
    return False

def partially_contains_sublist(main_list, list_of_lists):
    """
    Check if elements of the main list are partially contained in any sublist of the list of lists.

    Parameters
    ----------
    main_list : list
        The query list to check for partial containment.
    list_of_lists : list
        The subject list containing sublists to compare against the query list.

    Returns
    -------
    return : bool
        True if any element of the main list is partially contained in any sublist of the list of lists, False otherwise.
    """
    for sublist in list_of_lists:
        if any_match_lists(main_list, sublist):
            return True
    return False


def identify_string_in_dict(input_str, dictionary):
    """
    Identify the key in the dictionary where the input string is present.

    Parameters
    ----------
    input_str : str
        The string to find.
    dictionary : dict
        The dictionary where to find the string.

    Returns
    -------
    key : string or int
        The key of the entry where the string is present, or None if not found.
    """
    for key, value in dictionary.items():
        if input_str in value:
            return key
    return None

def has_element_of_type(input_list, target_type):
    """
    Check if any element in the input list matches the specified target type.

    Parameters
    ----------
    input_list : list
        The list to check for the presence of the target type.
    target_type : type
        The type to search for within the input list.

    Returns
    -------
    return : bool
        True if the target type is found in the input list, False otherwise.
    """
    for element in input_list:
        if isinstance(element, target_type):
            return True
    return False

def create_whitespace_string(input_string):
    """
    Create a string with the same length as the input string filled with whitespace.

    Parameters
    ----------
    input_string : str
        Input string to get the size from.

    Returns
    -------
    whitespace_string : str
        String containing the same number of white spaces as the length of the input string.
    """
    whitespace_string = " " * len(input_string)
    return whitespace_string

def remove_empty_dicts_recursive(nested_dict):
    """
    Recursively removes empty dictionary entries from a nested dictionary.

    Parameters
    ----------
    nested_dict : dict
        The nested dictionary.

    Returns
    -------
    nested_dict : dict
        The nested dictionary with empty dictionaries removed.
    """
    if isinstance(nested_dict, dict):
        # Iterate through a copy of the keys to avoid dictionary size change during iteration
        for key in list(nested_dict.keys()):
            nested_dict[key] = remove_empty_dicts_recursive(nested_dict[key])
            # Remove empty dictionaries (not considering numeric values like 0 as empty)
            if not nested_dict[key] and not isinstance(nested_dict[key], (bool, int, float)):
                del nested_dict[key]
    return nested_dict

def tsv_to_dict(file_path, skip_header = False, sep = '\t'):
    """
    Converts input TSV into a dict based on the desired separator.
    
    Parameters
    ----------
    file_path : str
        Path to the TSV file.
    skip_header : bool, optional
        If to ignore the header when importing.
    sep : str, optional
        String to seperate the values inside the TSV file.
    
    Returns
    -------
    data_dict : dict
        Returns the TSV file converted to the dict.
    """
    # Initialize an empty dictionary to store data
    data_dict = {}

    # Open the TSV file for reading
    with open(file_path, 'r') as f:
        # Read each line in the file
        for i, line in enumerate(f):
            if i == 0:
                continue
                
            # Split the line by tabs
            values = line.strip().split(sep)

            # Extract the key (first column entry)
            key = values[0]

            # Extract the rest of the values (excluding the first column entry)
            rest_values = values[1:]

            # Add the key-value pair to the dictionary
            data_dict[key] = rest_values

    return data_dict

def partially_contains_fragment_of_list(target_list, list_of_lists):
    """
    Check if the target_list is contained inside sublist even if it partially.
    e.g partially_contains_fragment_of_list(['a', 'b'], [['a', 'b', 'c'], ['d', 'e']])
    returns True.
    
    Parameters
    ----------
    target_list : list
        List to find inside the list_of_lists.
    list_of_lists : list
        The nested list.
    Returns
    -------
    returns : bool
        True if contains False if not.
    """
    for sub in list_of_lists:
        if any(sub[i:i+len(target_list)] == target_list for i in range(len(sub)-len(target_list)+1)):
            return True
    return False