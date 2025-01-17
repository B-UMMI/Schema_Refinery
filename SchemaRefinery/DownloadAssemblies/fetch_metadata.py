#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import time
import argparse
import concurrent.futures
from itertools import repeat
from typing import Any, Dict, List, Union
import xml.etree.ElementTree as ET

import pandas as pd
from tqdm import tqdm
from Bio import Entrez

try:
    from utils import (print_functions as pf)
except ModuleNotFoundError:
    from SchemaRefinery.utils import (print_functions as pf)

def find_internal_id(query_id: str) -> Union[str, int]:
    """
    To extract all metadata from the chosen Biosample it is needed to obtain
    internal id for all of the Biosamples.

    This function with an input of Biosample id fetches the internal id in NCBI
    for that Biosample.

    Parameters
    ----------
    query_id : str
        Biosample id

    Returns
    -------
    Union[str, int]
        Internal id used by NCBI entrez or "Failed" if an error occurs.
    """
    try:
        # Perform a search in the "biosample" database for the given query_id
        handle = Entrez.esearch(db="biosample", term=query_id)
        
        # Parse the XML response from the search
        tree = ET.parse(handle)
        root = tree.getroot()

        # Iterate over all "IdList" elements in the XML tree
        for id_list in root.iter("IdList"):
            # Iterate over all child elements of the current "IdList" element
            for in_id in id_list:
                # Check if the text content of the "Id" element is not None
                if in_id.text is not None:
                    # Convert the text content to an integer and return it as the internal ID
                    return int(in_id.text)
    except Exception:
        # Return "Failed" if an exception occurs during the process
        return "Failed"
    
    # Return "Failed" if no internal id is found
    return "Failed"

def download_query(internal_id: Union[str, int]) -> Union[ET.Element, str]:
    """
    This function with an input of an internal id for a Biosample, fetches
    the metadata related to that Biosample in xml object tree.

    Parameters
    ----------
    internal_id : Union[str, int]
        Internal id used by NCBI entrez

    Returns
    -------
    Union[ET.Element, str]
        XML tree object or "Failed" if an error occurs.
    """
    try:
        if internal_id != "Failed":
            # Fetch the metadata for the given internal_id from the "biosample" database
            handle = Entrez.efetch(db="biosample", id=str(internal_id))
            
            # Parse the XML response from the fetch
            tree = ET.parse(handle)
            root = tree.getroot()
            
            # Return the root element of the XML tree
            return root
        else:
            # Return "Failed" if the internal_id is "Failed"
            return "Failed"
    except Exception:
        # Return "Failed" if an exception occurs during the process
        return "Failed"

def get_metadata(xml_root: Union[ET.Element, str]) -> Dict[str, Any]:
    """
    This function based on input of an xml object tree, extracts relevant
    metadata inside a dictionary, where key is column and value is metadata.

    Parameters
    ----------
    xml_root : Union[ET.Element, str]
        XML tree object

    Returns
    -------
    Dict[str, Any]
        Dictionary that contains metadata
    """
    if isinstance(xml_root, ET.Element):
        # Initialize an empty dictionary to store metadata
        metadata_dict: Dict[str, Any] = {}
        
        # Update the dictionary with attributes from the root element
        metadata_dict.update(xml_root[0].attrib)
        
        # Find the "Organism" element in the XML tree
        organism_element = xml_root.find(".//Organism")
        if organism_element is not None:
            # Update the dictionary with attributes from the "Organism" element
            metadata_dict.update(organism_element.attrib)
        else:
            # Print a message if the "Organism" element is not found
            pf.print_message("Organism element not found.", "warning")

        # Iterate over all "Attributes" elements in the XML tree
        for attributes in xml_root.iter('Attributes'):
            # Iterate over all child elements of the current "Attributes" element
            for metadata in attributes:
                # Get the attribute name
                keyname: str = metadata.attrib['attribute_name']
                # Check if there is a harmonized name and use it as the key name
                if 'harmonized_name' in metadata.attrib:
                    keyname = metadata.attrib['harmonized_name']
                # Add the metadata to the dictionary
                metadata_dict[keyname] = metadata.text

        # Return the dictionary containing the metadata
        return metadata_dict
    else:
        # Return a dictionary with "accession" set to "Failed" if xml_root is not an ET.Element
        return {"accession": "Failed"}

def write_to_file(metadata_df: pd.DataFrame, output_directory: str) -> None:
    """
    Writes dataframe into TSV file format at output directory.

    Parameters
    ----------
    metadata_df : pd.DataFrame
        Pandas dataframe object containing metadata.
    output_directory : str
        Path to the output directory.

    Returns
    -------
    None
        Creates file in the output directory.
    """
    # Create the full path for the output file in the specified output directory
    output_file: str = os.path.join(output_directory, "metadata_biosamples.tsv")

    # Write the DataFrame to a TSV file
    # - mode='a' appends to the file if it exists, otherwise creates a new file
    # - header=not os.path.exists(output_file) writes the header only if the file does not already exist
    # - sep="\t" specifies that the file should be tab-separated
    # - index=False ensures that the DataFrame index is not written to the file
    metadata_df.to_csv(output_file, mode='a', header=not os.path.exists(output_file), sep="\t", index=False)

def multi_thread_run(query: str, retry: int) -> Dict[str, Any]:
    """
    Fetch metadata with multithreading.

    Parameters
    ----------
    query : str
        Biosample id
    retry : int
        Number of retries if fetching fails

    Returns
    -------
    Dict[str, Any]
        Dictionary that contains metadata
    """
    # Initialize the retry counter
    rtry: int = 0

    # Initialize an empty dictionary to store metadata
    metadata_dict: Dict[str, Any] = {}

    # Loop until the maximum number of retries is reached
    while rtry < retry:
        # Fetch the metadata by finding the internal ID, downloading the query, and extracting the metadata
        metadata_dict = get_metadata(download_query(find_internal_id(query)))

        # Check if the metadata dictionary contains only one item (indicating a failure)
        if len(metadata_dict) == 1:
            # Increment the retry counter
            rtry += 1
            # Wait for 1 second before retrying
            time.sleep(1)
        else:
            # Break the loop if the metadata dictionary contains more than one item (indicating success)
            break

    # Return the metadata dictionary
    return metadata_dict

def main(id_table_path: str, output_directory: str, email: str, threads: int, api_key: Union[str, None], retry: int) -> None:
    """
    Main function to handle fetching metadata based on provided arguments.

    Parameters
    ----------
    id_table_path : str
        Path to the input table containing Biosample ids.
    output_directory : str
        Path to the directory where the output files will be saved.
    email : str
        Email address for NCBI Entrez.
    threads : int
        Number of threads to use for fetching metadata.
    api_key : Union[str, None]
        API key to increase the number of requests.
    retry : int
        Number of retries if fetching fails.

    Returns
    -------
    None
        The function writes the output files to the specified directory.
    """
    Entrez.email = email

    # API key to increase number of requests
    if api_key is not None:
        Entrez.api_key = api_key

    # List where all the dictionaries containing metadata are stored
    metadata_list_dict: List[Dict[str, Any]] = []

    # Read input table to extract ids
    with open(id_table_path, 'r') as infile:
        queries: List[str] = infile.read().splitlines()

    failures: List[Dict[str, Any]] = []
    # Multithreading function
    with concurrent.futures.ThreadPoolExecutor(max_workers=threads) as executor:
        for res in list(tqdm(executor.map(multi_thread_run, queries, repeat(retry)), total=len(queries))):
            if 'Failed' in res:
                failures.append(res)
            metadata_list_dict.append(res)

    # Convert list of dictionaries to Pandas dataframe
    # Keys are column ids
    metadata_df: pd.DataFrame = pd.DataFrame(metadata_list_dict)

    # Organize metadata and associate with input id
    metadata_df.insert(0, "File", queries)

    # Write to TSV in output directory
    write_to_file(metadata_df, output_directory)


def parse_arguments():

    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument('-i', '--input_id_table', type=str,
                        required=True,
                        dest='id_table_path',
                        help='path to table with ID')

    parser.add_argument('-o', '--output-directory', type=str,
                        required=True,
                        dest='output_directory',
                        help='Path to the output directory.')

    parser.add_argument('-e', '--email', type=str,
                        required=True,
                        dest='email',
                        help='email for entrez queries')

    parser.add_argument('-t', '--threads', type=int,
                        required=False,
                        dest='threads',
                        default = 1,
                        help='number of threads, due to NCBI requests limitations recomended number of 2 without API key')

    parser.add_argument('-k', '--api_key', type=str,
                        required=False,
                        dest='api_key',
                        default=None,
                        help='API key to increase the mumber of requests')

    parser.add_argument('-r', '--retry', type=int,
                        required=False, dest='retry',
                        default=7,
                        help='Maximum number of retries when a '
                             'download fails.')

    args = parser.parse_args()

    return args


if __name__ == '__main__':

    args = parse_arguments()
    main(**vars(args))
