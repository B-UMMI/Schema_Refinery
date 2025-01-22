from collections import OrderedDict, Counter
import re
import itertools
import pickle
from typing import List, Tuple, Dict, Union, Any, Set


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


def decode_CDS_sequences_ids(path_to_file: str) -> Dict[str, List[float]]:
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


def identify_string_in_dict_get_key(input_str: str, dictionary: Dict[Union[str, int], Any]) -> Union[str, int, None]:
    """
    Identify the key in the dictionary where the input string is present.

    Parameters
    ----------
    input_str : str
        The string to find.
    dictionary : Dict[Union[str, int], Any]
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


def get_shared_elements(dict_: Dict[str, Set[Any]]) -> List[Any]:
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
    # Flatten all sets into a single list
    all_elements = [elem for subset in dict_.values() for elem in subset]
    
    # Count the occurrences of each element
    element_counts = Counter(all_elements)
    
    # Identify elements that appear in at least two sets
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
