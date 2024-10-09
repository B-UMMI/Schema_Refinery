from collections import OrderedDict, Counter
import re
import itertools
import pickle
from typing import List, Tuple, Dict, Union, Any, Set

def join_list(lst: List[Any], delimiter: str) -> str:
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
    str
        A single string with all elements in the input list joined by the specified delimiter.
    """
    return delimiter.join(map(str, lst))

def flatten_list(list_to_flatten: List[List[Any]]) -> List[Any]:
    """
    Flatten one level of a nested list.

    Parameters
    ----------
    list_to_flatten : list
        Nested list to flatten.

    Returns
    -------
    list
        Input list flattened by one level.
    """
    return list(itertools.chain(*list_to_flatten))

def isListEmpty(input_list: List[Any]) -> bool:
    """
    Check if a nested list is empty.
    
    Parameters
    ----------
    input_list : list
        Input list to verify

    Returns
    -------
    bool
        If list is empty or not.
    """
    if isinstance(input_list, list):
        return all(map(isListEmpty, input_list)) if isinstance(input_list, list) else False

def divide_list_into_n_chunks(list_to_divide: List[Any], n: int) -> List[List[Any]]:
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
    list
        A list of sublists created by dividing the input list.
    """
    sublists = []
    list_length = len(list_to_divide)
    sublist_size, remainder = divmod(list_length, n)
    
    start_idx = 0
    for i in range(n):
        sublist_length = sublist_size + (1 if i < remainder else 0)
        end_idx = start_idx + sublist_length
        sublists.append(list_to_divide[start_idx:end_idx])
        start_idx = end_idx

    return [sublist for sublist in sublists if sublist]

def get_max_min_values(input_list: List[Union[int, float]]) -> Tuple[Union[int, float], Union[int, float]]:
    """
    From an input list, return the maximum and minimum integer or float values.

    Parameters
    ----------
    input_list : list[int] or list[float]
        List containing integer or float values.

    Returns
    -------
    tuple
        A tuple containing the largest and smallest values in the input list.
    """
    return max(input_list), min(input_list)

def decompress_number(text: str, index: int) -> Tuple[int, int]:
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
    tuple
        Index to start decoding the next number and the number decoded
        in the current function call.
    """
    number = 0
    bitwise_shift = 0

    while True:
        n = (ord(text[index]) - 63)
        index += 1
        if n >= 0x20:
            n -= 0x20
            number = number | (n << bitwise_shift)
            bitwise_shift += 5
        else:
            break

    number = number | (n << bitwise_shift)
    return index, (~number >> 1) if (number & 1) != 0 else (number >> 1)

def polyline_decoding(text: str) -> List[float]:
    """
    Decode a list of integers compressed with polyline encoding.

    Parameters
    ----------
    text : str
        String representing a compressed list of numbers.

    Returns
    -------
    list
        List with the decoded numbers.
    """
    number_list = []
    index = last_num = 0
    precision = ord(text[index]) - 63
    index += 1
    while index < len(text):
        index, difference = decompress_number(text, index)
        last_num += difference
        number_list.append(last_num)

    return [round(item * (10 ** (-precision)), precision) for item in number_list]

def decode_CDS_sequences_ids(path_to_file: str) -> Dict[str, List[int]]:
    """
    Read a dictionary contained in a pickle file and decode its values based on polyline.

    Parameters
    ----------
    path_to_file : str
        Path to the pickle file.

    Returns
    -------
    dict
        Contains hashed CDS as keys, and number id of the genome where that CDS is present.
    """
    with open(path_to_file, "rb") as infile:
        hash_table = pickle.load(infile)

    return {key: polyline_decoding(value) for key, value in hash_table.items()}

def get_unique_sublists(list_of_lists: List[List[Any]]) -> List[List[Any]]:
    """
    Identify unique sublists within a list of lists.

    Parameters
    ----------
    list_of_lists : list
        The list containing various sublists.

    Returns
    -------
    list
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

def all_match_lists(list1: List[Any], list2: List[Any]) -> bool:
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
    bool
        True if all elements of list1 are found in list2, False otherwise.
    """
    return all(elem in list2 for elem in list1)

def any_match_lists(list1: List[Any], list2: List[Any]) -> bool:
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
    bool
        True if any element of list1 is found in list2, False otherwise.
    """
    return any(elem in list2 for elem in list1)

def contains_sublist(main_list: List[Any], list_of_lists: List[List[Any]]) -> bool:
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
    bool
        True if all elements of the main list are fully contained in any sublist of the list of lists, False otherwise.
    """
    for sublist in list_of_lists:
        if all_match_lists(main_list, sublist):
            return True
    return False

def partially_contains_sublist(main_list: List[Any], list_of_lists: List[List[Any]]) -> bool:
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
    bool
        True if any element of the main list is partially contained in any sublist of the list of lists, False otherwise.
    """
    for sublist in list_of_lists:
        if any_match_lists(main_list, sublist):
            return True
    return False

def identify_string_in_dict_get_key(input_str: str, dictionary: Dict[Union[str, int], Any]) -> Union[str, int, None]:
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
    str or int or None
        The key of the entry where the string is present, or None if not found.
    """
    for key, value in dictionary.items():
        if input_str in value:
            return key
    return None

def identify_string_in_dict_get_value(search_key: str, search_dict: Dict[str, Any]) -> Any:
    """
    Identify the value in the dictionary where the search key is present.

    Parameters
    ----------
    search_key : str
        The key to find.
    search_dict : dict
        The dictionary where to find the key.

    Returns
    -------
    Any
        The value of the entry where the key is present, or None if not found.
    """
    for key, value in search_dict.items():
        if search_key in value:
            return value
    return None

def identify_string_in_dict_lists_regex(target_value: str, dict_of_lists: Dict[Union[str, int], List[List[str]]], regex: bool = False) -> Union[str, int, bool]:
    """
    Identifies if a string is present in any list inside a dictionary.

    Parameters
    ----------
    target_value : str
        The value to find.
    dict_of_lists : dict
        A dictionary where the values are lists of lists.
    regex : bool, optional
        A regex pattern to search for in the lists.

    Returns
    -------
    int or str or bool
        The key of the entry where the string is present, or False if not found.
    """
    for key, lists in dict_of_lists.items():
        for list_ in lists:
            if target_value in [remove_by_regex(l, regex) for l in list_] if regex else list_:
                return key
    return False

def has_element_of_type(input_list: List[Any], target_type: type) -> bool:
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
    bool
        True if the target type is found in the input list, False otherwise.
    """
    for element in input_list:
        if isinstance(element, target_type):
            return True
    return False

def create_whitespace_string(input_string: str) -> str:
    """
    Create a string with the same length as the input string filled with whitespace.

    Parameters
    ----------
    input_string : str
        Input string to get the size from.

    Returns
    -------
    str
        String containing the same number of white spaces as the length of the input string.
    """
    return " " * len(input_string)

def remove_empty_dicts_recursive(nested_dict: Dict[Any, Any]) -> Dict[Any, Any]:
    """
    Recursively removes empty dictionary entries from a nested dictionary.

    Parameters
    ----------
    nested_dict : dict
        The nested dictionary.

    Returns
    -------
    dict
        The nested dictionary with empty dictionaries removed.
    """
    if isinstance(nested_dict, dict):
        for key in list(nested_dict.keys()):
            nested_dict[key] = remove_empty_dicts_recursive(nested_dict[key])
            if not nested_dict[key] and not isinstance(nested_dict[key], (bool, int, float)):
                del nested_dict[key]
    return nested_dict

def tsv_to_dict(file_path: str, skip_header: bool = False, sep: str = '\t') -> Dict[str, List[str]]:
    """
    Converts input TSV into a dict based on the desired separator.
    
    Parameters
    ----------
    file_path : str
        Path to the TSV file.
    skip_header : bool, optional
        If to ignore the header when importing.
    sep : str, optional
        String to separate the values inside the TSV file.
    
    Returns
    -------
    dict
        Returns the TSV file converted to the dict.
    """
    data_dict = {}

    with open(file_path, 'r') as f:
        for i, line in enumerate(f):
            if i == 0 and skip_header:
                continue
                
            values = line.strip().split(sep)
            key = values[0]
            rest_values = values[1:]
            data_dict[key] = rest_values

    return data_dict

def partially_contains_fragment_of_list(target_list: List[Any], list_of_lists: List[List[Any]]) -> bool:
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
    bool
        True if contains False if not.
    """
    for sub in list_of_lists:
        if any(sub[i:i+len(target_list)] == target_list for i in range(len(sub)-len(target_list)+1)):
            return True
    return False

def remove_by_regex(string: str, pattern: str) -> str:
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
    str
        The string with all occurrences of the pattern removed.
    """
    return re.sub(pattern, '', string)

def replace_by_regex(string: str, pattern: str, replacement: str) -> str:
    """
    Replace all occurrences of a regex pattern within a string with a specified replacement.

    Parameters
    ----------
    string : str
        The string to search and replace occurrences in.
    pattern : str
        The regex pattern to search for within `string`.
    replacement : str
        The string to replace each match of `pattern` in `string` with.

    Returns
    -------
    str
        A new string with all matches of `pattern` replaced by `replacement`.
    """
    return re.sub(pattern, replacement, string)

def regex_present(regex_list: List[str], string: str) -> bool:
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
    bool
        True if any regex is found in the string, False otherwise.
    """
    return any(re.search(regex, string) for regex in regex_list)

def search_string_by_regex(pattern: str, string: str) -> Union[str, None]:
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
    str or None
        The match object if the pattern is found, original string otherwise.
    """
    match = re.search(pattern, string)
    return match.group(1) if match else string

def add_strings_to_subsets(my_list: List[Set[str]], my_strings: List[str]) -> bool:
    """
    Clustering algorithm that finds a string in a list of strings in
    a list of sets and adds the whole list to the set if any string of 
    that list is inside the set.

    Parameters
    ----------
    my_list : list of sets
        The list of sets to add strings to.
    my_strings : list of str
        The list of strings to add to the sets.

    Returns
    -------
    bool
        True if any string was added to a set, False otherwise.
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

def find_index(input_list: List[str], target_string: str) -> Union[int, None]:
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
    int or None
        The index of the string in the list, or None if the string is not found.
    """
    try:
        return input_list.index(target_string)
    except ValueError:
        return None

def find_sublist_index(input_list_of_lists: List[List[Any]], target_value: Any) -> Union[int, None]:
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
    int or None
        The index of the sublist containing the element, or None if the string is not found in any sublist.
    """
    try:
        return input_list_of_lists.index(target_value)
    except ValueError:
        return None

def try_convert_to_type(value: Any, target_type: type) -> Any:
    """
    Attempts to convert a given value to a specified type.

    Parameters
    ----------
    value : Any
        The value to be converted.
    target_type : type
        The type to which `value` should be converted.

    Returns
    -------
    Any
        The converted value if the conversion is successful; otherwise, the original value.
    """
    try:
        return target_type(value)
    except (ValueError, TypeError):
        return value

def repeat_strings_in_a_list(string: str, times: int) -> List[str]:
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
    """
    return [string for _ in range(times)]

def sort_subdict_by_tuple(dict_: Dict[str, Dict[str, Any]], order: Tuple[str, ...]) -> Dict[str, OrderedDict]:
    """
    Sorts the sub-dictionaries of a given dictionary based on a specified order tuple.

    Parameters
    ----------
    dict_ : dict
        The input dictionary containing sub-dictionaries as values.
    order : tuple
        A tuple specifying the desired order of keys in the sorted sub-dictionaries.

    Returns
    -------
    dict
        A new dictionary with each sub-dictionary sorted according to the specified order.
    """
    sorted_data = {}
    for key, subdict in dict_.items():
        sorted_subdict = OrderedDict(sorted(subdict.items(), key=lambda item: order.index(item[0]) if item[0] in order else len(order)))
        sorted_data[key] = sorted_subdict
    return sorted_data

def check_if_all_elements_are_duplicates(input_list: List[Any]) -> bool:
    """
    Check if all elements in the list are duplicates.

    Parameters
    ----------
    input_list : list
        The list to check for duplicate elements.

    Returns
    -------
    bool
        True if every element in the list occurs more than once, False otherwise.
        Returns False if the list is empty.
    """
    element_counts = {}
    for element in input_list:
        if element in element_counts:
            element_counts[element] += 1
        else:
            element_counts[element] = 1
    
    for count in element_counts.values():
        if count == 1:
            return False
    return True if element_counts else False

def check_if_all_sets_are_same(sets_list: List[Set[Any]]) -> bool:
    """
    Checks if all sets within a list are identical.

    Parameters
    ----------
    sets_list : list
        A list of sets to be checked for identity.

    Returns
    -------
    bool
        True if all sets in the list are identical, False otherwise.
    """
    if len(sets_list) <= 1:
        return True
    
    reference_set = sets_list[0]
    
    for s in sets_list[1:]:
        if s != reference_set:
            return False
    
    return True

def get_duplicates(input_list: List[Any]) -> List[Any]:
    """
    Identify duplicate elements in a list.

    Parameters
    ----------
    input_list : list
        The list to check for duplicate elements.

    Returns
    -------
    list
        A list containing the duplicate elements found in the input list.
    """
    element_counts = Counter(input_list)
    duplicates = [element for element, count in element_counts.items() if count > 1]
    return duplicates

def join_list(lst: List[str], delimiter: str) -> str:
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
    str
        A single string with all elements in the input list joined by the character chosen as link.
    """
    return delimiter.join(lst)

def get_common_elements_in_lists(list_of_lists: List[List[Any]]) -> List[Any]:
    """
    Finds common elements between various lists.

    Parameters
    ----------
    list_of_lists : list
        Contains a list of lists.

    Returns
    -------
    list
        Returns a list that contains the intersection of all elements inside the list of lists.
    """
    intersection_set = None
    for lst in list_of_lists:
        if intersection_set is None:
            intersection_set = set(lst)
        else:
            intersection_set.intersection_update(lst)

    return list(intersection_set)

def get_shared_elements(dict_: Dict[str, List[Any]]) -> List[Any]:
    """
    Identify elements that appear in at least two lists within a dictionary.

    Parameters
    ----------
    dict_ : dict
        A dictionary where the values are lists of elements.

    Returns
    -------
    list
        A list containing elements that appear in at least two lists within the dictionary.
    """
    all_elements = [elem for sublist in dict_.values() for elem in sublist]
    element_counts = Counter(all_elements)
    shared_elements = [elem for elem, count in element_counts.items() if count >= 2]
    return shared_elements

def convert_set_elements_to_strings(input_set: Set[Any]) -> Set[str]:
    """
    Convert all elements in the set to strings.

    Parameters
    ----------
    input_set : set
        The set of elements to be converted to strings.

    Returns
    -------
    set
        A set with all elements converted to strings.
    """
    return {str(element) for element in input_set}