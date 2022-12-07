#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 28 13:28:18 2022

@author: Mykyta Forofontov
"""

import os
import pandas as pd
import argparse
import xml.etree.ElementTree as ET
import concurrent.futures
import time
from Bio import Entrez
from itertools import repeat

def read_table(id_table_path):
    """
    Reads TSV file containing Biosample id and transforms into pandas dataframe.
    
    input: TSV file path
    output: pandas dataframe object
    """    
    with open(id_table_path) as id_table:
        
        return pd.read_csv(id_table, sep='\t',low_memory=False,header=None)
    
def find_internal_id(query_id):
    """
    To extract all metadata from the chosen Biosample it is needed to obtain
    internal id for all of the Biosamples.
    
    This function with an input of Biosample fetches the internal id
    for that Biosample.
    
    input: Biosample id
    output: internal id
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
    
    input: internal id
    output: xml tree object
    """
    
    try:
        if not "Failed" in internal_id:
                
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
    
    input: xml tree object
    output: dictionary
    """
    
    if not "Failed" in xml_root:
        
        metadata_dict = dict()
        
        metadata_dict.update(xml_root[0].attrib)
            
        metadata_dict.update(xml_root.find(".//Organism").attrib)
            
        metadata_dict["SRA"] = xml_root.findall(".//Id")[1].text

                          
        for attributes in xml_root.iter('Attributes'):
                
            for metadata in attributes:
                keyname = metadata.attrib['attribute_name']
                if 'harmonized_name' in metadata.attrib:
                    keyname = metadata.attrib['harmonized_name']
                metadata_dict[keyname] = metadata.text

        return metadata_dict
    
    else:
        return {"accession": "Failed"}
    
def write_to_file(metadata_df,output_directory):
    """
    Writes dataframe into TSV file format at output directory.
    
    input: pandas dataframe object
    output: TSV file inside output path
    """    
    metadata_df.to_csv(os.path.join(output_directory, "metadata.tsv"),
                       mode='a', header=not os.path.exists(os.path.join(
                           output_directory, "metadata.tsv")),sep="\t",index=False)
    
def multi_thread_run(query,retry):
    """
    Function that enables main function to run multithreading for faster metadata
    fetching
    
    input: Biosample id
    output: dictionary containg metadata for the input Biosample
    """
    
    print("\rDownloading {}".format(query),end=' ')

    rtry = 0
    
    while rtry < retry:
        
        metadata_dict = get_metadata(download_query(find_internal_id(query)))
        
        if len(metadata_dict) == 1:
            rtry += 1
            time.sleep(1)
        else:
            break
    
    return metadata_dict
            
def main(id_table_path,output_directory,email,threads,api_key,retry):
    
    #Create directory if absent in output path
    if not os.path.exists(output_directory):
        os.mkdir(output_directory)    
    
    Entrez.email = email
    
    #API key to increase number of requests
    if api_key != None:
        Entrez.api_key = api_key
    
    #list where all the dictionaries containg metadata are stored
    metadata_list_dict = list()
    
    #Read input table to extract ids
    table = read_table(id_table_path)
    queries = [i for sl in table.values.tolist() for i in sl]
    
    failures = list()
    
    #multithreading function
    with concurrent.futures.ThreadPoolExecutor(max_workers=threads) as executor:

        for res in executor.map(multi_thread_run, queries, repeat(retry)):
            if 'Failed' in res:
                failures.append(res)
            metadata_list_dict.append(res)
            
    """
    Convert the list of dictionaries into pandas dataframe where keys are
    columns and values are values for the following column.
    When dictionary does not have a key and value, the value in the column of
    pandas dataframe if that column already exists is kept blank.
    """
    
    metadata_df = pd.DataFrame(metadata_list_dict)
    
    #Organise metadata and associate with input id
    metadata_df.insert(0,"accession",metadata_df.pop("accession"))
    metadata_df.insert(0,"File",queries)
    
    #Write to TSV in output directory
    write_to_file(metadata_df,output_directory)

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
