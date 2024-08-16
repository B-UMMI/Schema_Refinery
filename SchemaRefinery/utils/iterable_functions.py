from collections import OrderedDict, Counter
import re
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


def identify_string_in_dict_get_key(input_str, dictionary):
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

def identify_string_in_dict_get_value(search_key, search_dict):
    for key, value in search_dict.items():
        if search_key in value:
            return value
    return None

def identify_string_in_dict_lists_regex(target_value, dict_of_lists, regex=False):
    """
    Identifies if a string is present in any list inside a dictionary.

    Parameters
    ----------
    target_string : str
        The value to find.
    dict_of_lists : dict
        A dictionary where the values are lists of lists.
    regex : str, optional
        A regex pattern to search for in the lists.

    Returns
    -------
    key : int or str
        The key of the entry where the string is present, or False if not found.
    """
    for key, lists in dict_of_lists.items():
        for list_ in lists:
            if target_value in [remove_by_regex(l, regex) for l in list_] if regex else list_:
                return key
    return False

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

def remove_by_regex(string, pattern):
    """
    Remove all occurrences of a pattern from a string.

    Parameters
    ----------
    string : str
        The string to remove the pattern from.
    pattern : str
        The regex pattern to remove from the string.

    Returns
    -------
    return : str
        The string with all occurrences of the pattern removed.
    """
    return re.sub(pattern, '', string)

def replace_by_regex(string, pattern, replacement):
    """
    Replace all occurrences of a regex pattern within a string with a specified replacement.

    This function uses the `re.sub()` method from Python's built-in `re` (regular expressions) module to find all
    occurrences of `pattern` in `string` and replace them with `replacement`. The function returns a new string
    with the modifications applied.

    Parameters
    ----------
    string : str
        The string to search and replace occurrences in.
    pattern : str
        The regex pattern to search for within `string`. This pattern can match characters, numbers, symbols,
        or more complex regex features.
    replacement : str
        The string to replace each match of `pattern` in `string` with.

    Returns
    -------
    str
        A new string with all matches of `pattern` replaced by `replacement`.

    Examples
    --------
    >>> text = "Hello 123, meet 456."
    >>> pattern = r"\d+"
    >>> replacement = "number"
    >>> replace_by_regex(text, pattern, replacement)
    'Hello number, meet number.'
    """
    return re.sub(pattern, replacement, string)

def regex_present(regex_list, string):
    """
    Check if any regex in a list is found in a string.

    Parameters
    ----------
    regex_list : list of str
        The list of regexes to search for in the string.
    string : str
        The string to search in.

    Returns
    -------
    return : bool
        True if any regex is found in the string, False otherwise.
    """
    return any(re.search(regex, string) for regex in regex_list)

def search_string_by_regex(pattern, string):
    """
    Searches for a regex pattern in a string.

    Parameters
    ----------
    pattern : str
        The regex pattern to search for.
    string : str
        The string to search in.

    Returns
    -------
    return : str
        The match object if the pattern is found, original string otherwise.
    """
    match = re.search(pattern, string)
    return match.group(1) if match else string

def add_strings_to_subsets(my_list, my_strings):
    """
    Clustering algorithm that finds a string in a list of strings in
    a list of sets and adds th whole list to the set if any string of 
    that list is inside the set
    """
    found = False
    for my_string in my_strings:
        if found:
            break
        for sublist in my_list:
            if my_string in sublist:
                sublist.update(my_strings)
                found = True
                break
    return found
#Unused
def find_index(input_list, target_string):
    """
    Find the index of a string in a list.

    Parameters
    ----------
    input_list : list
        The list to search.
    target_string : str
        The string to find.

    Returns
    -------
    index : int
        The index of the string in the list, or None if the string is not found.
    """
    try:
        return input_list.index(target_string)
    except ValueError:
        return None

def find_sublist_index(input_list_of_lists, target_value):
    """
    Finds the index of the sublist that contains the target value within a list of lists.

    Parameters
    ----------
    input_list_of_lists : list of list
        The list of lists to search.
    target_value : any
        The value to find.

    Returns
    -------
    return : int or None
        The index of the sublist containing the element, or None if the string is not found in any sublist.
    """
    try:
        return input_list_of_lists.index(target_value)
    except ValueError:
        return None

def try_convert_to_type(value, target_type):
    """
    Attempts to convert a given value to a specified type.

    This function tries to convert the input `value` to the `target_type`. If the conversion is successful, it returns the converted value. If the conversion fails due to a `ValueError` or `TypeError`, it returns the original value instead.

    Parameters
    ----------
    value : Any
        The value to be converted.
    target_type : type
        The type to which `value` should be converted. This should be a type like `int`, `float`, `str`, etc.

    Returns
    -------
    return : Any
        The converted value if the conversion is successful; otherwise, the original value.

    Notes
    -----
    - This function is useful for safely attempting type conversions without the risk of raising exceptions.
    - It can be used in situations where the type of input data is uncertain or varies.
    """
    try:
        return target_type(value)
    except (ValueError, TypeError):
        return value
#Unused
def repeat_strings_in_a_list(string, times):
    """
    Creates a list where a given character is repeated a specified number of times.

    Parameters
    ----------
    string : str
        The character to be repeated.
    times : int
        The number of times the character should be repeated.

    Returns
    -------
    list of str
        A list containing the character repeated 'times' times.

    Examples
    --------
    >>> repeat_strings_in_a_list('a', 3)
    ['a', 'a', 'a']
    """
    return [string for i in range(times)]

def sort_subdict_by_tuple(dict, order):
    """
    Sorts the sub-dictionaries of a given dictionary based on a specified order tuple.

    Parameters
    ----------
    dict : dict
        The input dictionary containing sub-dictionaries as values.
    order : tuple
        A tuple specifying the desired order of keys in the sorted sub-dictionaries.

    Returns
    -------
    sorted_data : dict
        A new dictionary with each sub-dictionary sorted according to the specified order.

    Notes
    -----
    -This function iterates through each key-value pair in the input dictionary. Each value, 
    which should be a dictionary itself (sub-dictionary), is sorted based on the order of keys 
    specified in the 'order' tuple. If a key in the sub-dictionary does not exist in the 'order' tuple, 
    it is placed at the end of the sorted sub-dictionary. The sorting is stable, meaning that 
    the original order of keys (for those not in the 'order' tuple) is preserved.

    Examples
    --------
    >>> data = {'cds1|cds2': {'1a': 30, '3b': 40}, 'cds3|cds4': {'3b': 20, '1a': 30}}
    >>> order = ('1a', '1b', '2a', '3a', '2b', '1c', '3b', '4a', '4b', '4c', '5')
    >>> sorted_data = sort_subdict_by_tuple(data, order)
    >>> sorted_data
    {'cds1|cds2': OrderedDict([('1a', 30), ('3b', 40)]), 'cds3|cds4': OrderedDict([('1a', 30), ('3b', 20)])}
    """
    sorted_data = {}
    for key, subdict in dict.items():
        # Sorting the sub-dictionary by the index of its keys in the order tuple
        sorted_subdict = OrderedDict(sorted(subdict.items(), key=lambda item: order.index(item[0]) if item[0] in order else len(order)))
        sorted_data[key] = sorted_subdict
    return sorted_data

def check_if_all_elements_are_duplicates(input_list):
    # Count occurrences of each element
    element_counts = {}
    for element in input_list:
        if element in element_counts:
            element_counts[element] += 1
        else:
            element_counts[element] = 1
    
    # Check if every element occurs more than once
    for count in element_counts.values():
        if count == 1:
            return False
    return True if element_counts else False

def check_if_all_sets_are_same(sets_list):
    """
    Checks if all sets within a list are identical.

    This function evaluates whether all sets in a given list are exactly the same. It first checks if the list is empty or contains only one set, in which case it returns True, as there are no differing sets to compare. Then, it uses the first set in the list as a reference to compare against all other sets in the list. If any set differs from the first set, the function returns False. Otherwise, if all sets are identical to the first set, it returns True.

    Parameters
    ----------
    sets_list : list
        A list of sets to be checked for identity.

    Returns
    -------
    bool
        True if all sets in the list are identical, False otherwise.

    Examples
    --------
    >>> sets_list = [{1, 2, 3}, {1, 2, 3}, {1, 2, 3}]
    >>> check_if_all_sets_are_same(sets_list)
    True

    >>> sets_list2 = [{1, 2, 3}, {4, 5, 6}, {1, 2, 3}]
    >>> check_if_all_sets_are_same(sets_list2)
    False
    """
    # Check if the list is empty or has only one set
    if len(sets_list) <= 1:
        return True
    
    # Use the first set as a reference
    reference_set = sets_list[0]
    
    # Compare each set with the reference set
    for s in sets_list[1:]:
        if s != reference_set:
            return False
    
    return True

def get_duplicates(input_list):
    # Count occurrences of each element
    element_counts = Counter(input_list)
    # Select elements that appear more than once
    duplicates = [element for element, count in element_counts.items() if count > 1]
    return duplicates