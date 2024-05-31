import pandas as pd

def dict_to_df(dictionary):
    """
    Convert a dictionary to a pandas DataFrame.

    Parameters
    ----------
    dictionary : dict
        The dictionary to convert. Keys are used as column headers.

    Returns
    -------
    return : pandas.DataFrame
        The resulting DataFrame, where each key-value pair in the dictionary corresponds to a column in the DataFrame.
    """
    return pd.DataFrame.from_dict(dictionary)