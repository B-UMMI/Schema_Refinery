#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
"""

import os
import time
import argparse
import concurrent.futures
from itertools import repeat
import xml.etree.ElementTree as ET

import pandas as pd
from tqdm import tqdm
from Bio import Entrez


def find_internal_id(query_id):
    """
    To extract all metadata from the chosen Biosample it is needed to obtain
    internal id for all of the Biosamples.

    This function with an input of Biosample id fetches the internal id in NCBI
    for that Biosample.

    Parameter
    ---------
    query_id : str
        Biosample id

    Returns
    -------
    internal_id : int
        internal id used by NCBI entrez
    """
    try:
        handle = Entrez.esearch(db="biosample", term=query_id)
        tree = ET.parse(handle)
        root = tree.getroot()

        for id_list in root.iter("IdList"):

            for in_id in id_list:

                return in_id.text
    except:
        return "Failed"


def download_query(internal_id):
    """
    This function with and input of an internal id for a Biosample, fetches
    the metadata related to that Biosample in xml object tree.

    Parameter
    ---------
    internal_id : int
        internal id used by NCBI entrez

    Returns
    -------
    root : object
        xml tree object
    internal_id : str
        internal id of the sample
    "Failed" : str
        if errors ocurrs return string "Failed"
    """
    try:
        if "Failed" not in internal_id:
            handle = Entrez.efetch(db="biosample", id=internal_id)
            tree = ET.parse(handle)
            root = tree.getroot()

            return root
        else:
            return internal_id
    except:
        return "Failed"


def get_metadata(xml_root):
    """
    This function based on input of an xml object tree, extracts relevant
    metadata inside a dictionary, where key is column and value is metadata

    Parameter
    ---------

    xml_root : xml tree object

    Returns
    -------
    return : dict
        that contains metadata
    """
    if not "Failed" in xml_root:
        metadata_dict = {}
        metadata_dict.update(xml_root[0].attrib)
        metadata_dict.update(xml_root.find(".//Organism").attrib)

        for attributes in xml_root.iter('Attributes'):

            for metadata in attributes:
                keyname = metadata.attrib['attribute_name']
                if 'harmonized_name' in metadata.attrib:
                    keyname = metadata.attrib['harmonized_name']
                metadata_dict[keyname] = metadata.text

        return metadata_dict

    else:
        return {"accession": "Failed"}


def write_to_file(metadata_df, output_directory):
    """
    Writes dataframe into TSV file format at output directory.

    Parameter
    ---------

    metadata_df : pandas dataframe object
    output_directory: str

    Returns
    -------
        Creates file in the output directory
    """
    metadata_df.to_csv(os.path.join(output_directory, "metadata.tsv"),
                       mode='a', header=not os.path.exists(os.path.join(
                           output_directory, "metadata.tsv")),sep="\t",index=False)


def multi_thread_run(query, retry):
    """Fetch metadata with multithreading.

    Parameter
    ---------
    query : str
        Biosample id
    retry : int
        number of retries if fetching fails

    Returns
    -------
    metadata_dict : dict
        dict that contains metadata
    """
    rtry = 0
    while rtry < retry:
        metadata_dict = get_metadata(download_query(find_internal_id(query)))

        if len(metadata_dict) == 1:
            rtry += 1
            time.sleep(1)
        else:
            break

    return metadata_dict


def main(id_table_path, output_directory, email, threads, api_key, retry):
    Entrez.email = email

    # API key to increase number of requests
    if api_key is not None:
        Entrez.api_key = api_key

    # list where all the dictionaries containg metadata are stored
    metadata_list_dict = []

    # Read input table to extract ids
    with open(id_table_path, 'r') as infile:
        queries = infile.read().splitlines()

    failures = []
    # multithreading function
    with concurrent.futures.ThreadPoolExecutor(max_workers=threads) as executor:
        for res in list(tqdm(executor.map(multi_thread_run, queries, repeat(retry)),total=len(queries))):
            if 'Failed' in res:
                failures.append(res)
            metadata_list_dict.append(res)

    # convert list of dctionaries to Pandas dataframe
    # keys are column ids
    metadata_df = pd.DataFrame(metadata_list_dict)

    # Organise metadata and associate with input id
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
