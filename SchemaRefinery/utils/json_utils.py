from typing import Any, Dict, List

def merge_keys(dicts: List[Dict[str, Any]]) -> Dict[str, Any]:
    """
    Recursively merges keys from a list of dictionaries into a single dictionary.

    This function takes a list of dictionaries and merges all the keys, including nested keys,
    into a single dictionary. The keys are represented in a dot-separated format to indicate
    the hierarchy.

    Parameters
    ----------
    dicts : List[Dict[str, Any]]
        A list of dictionaries to merge keys from.

    Returns
    -------
    Dict[str, Any]
        A dictionary containing all the merged keys from the input dictionaries.
    """
    merged: Dict[str, Any] = {}

    def merge(d: Dict[str, Any], parent_key: str = '') -> None:
        """
        Recursively merges keys from a dictionary into the merged dictionary.

        Parameters
        ----------
        d : Dict[str, Any]
            The dictionary to merge keys from.
        parent_key : str, optional
            The parent key to use for dot-separated key representation (default is '').
        """
        for key, value in d.items():
            # Create a dot-separated key if parent_key is provided
            full_key = f"{parent_key}.{key}" if parent_key else key
            if isinstance(value, dict):
                # Recursively merge keys from nested dictionaries
                merge(value, full_key)
            elif isinstance(value, list):
                # Handle lists of dictionaries
                if value and isinstance(value[0], dict):
                    for item in value:
                        merge(item, full_key)
                else:
                    # Handle lists of non-dictionary items
                    merged[full_key] = None
            else:
                # Add the key to the merged dictionary
                merged[full_key] = None

    # Merge keys from each dictionary in the list
    for d in dicts:
        merge(d)

    return merged
