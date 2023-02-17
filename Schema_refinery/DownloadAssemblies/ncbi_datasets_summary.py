#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Dec 27 13:05:38 2022

AUTHOR

    Mykyta Forofontov
    github: @MForofontov
"""

import subprocess
import json
import csv
import sys
import os

def verify_assembly(metadata_assembly,size_threshold,max_contig_number,
                    genome_size,verify_status):
    """
    This function verifies assemblies by certain inputa criteria.

    Parameters
    ----------
    metadata_assembly : json object (dict)
        For a single assembly.
    size_threshold : float
        (0 >= x >= 1).
    max_contig_number: int
        (>0).
    genome_size: int
        (>0).
    verify_status: bool

    Returns
    -------
    Boolean value (in order to see if passed or failed)
    """
    assembly_stats = metadata_assembly['assembly_stats']
    assembly_info = metadata_assembly['assembly_info']

    if genome_size is not None and size_threshold is not None:

        bot_limit = genome_size - (genome_size*size_threshold)
        top_limit = genome_size + (genome_size*size_threshold)

        if int(assembly_stats['total_sequence_length']) >= top_limit:

            return False

        if int(assembly_stats['total_sequence_length']) <= bot_limit:

            return False

    if max_contig_number is not None:

        if assembly_stats['number_of_contigs'] > max_contig_number:

            return False

    if verify_status is True or verify_status is None:
        if assembly_info['assembly_status'] == 'suppressed':
            return False

    return True


def metadata_fetcher_ids(input_ids,id_list_path,assembly_level,reference,
                         exclude_atypical,api_key):

    """
    This function based on an input id fetches json object (dict) for all
    assemblies.

    Parameters
    ----------
    input_id : list of str
        starts with GCF_ or GCA_ .
    assembly_level : str
        containing assembly levels separated by ",".
    reference: Bool
    api_key: str

    Returns
    -------
    found_metadata: bool
    failed_list: list
        containing ids that failed initial processing.
    metadata_assembly: json object (dict)
        for all assemblies.
    """

    arguments = ['datasets','summary','genome','accession', '--inputfile', id_list_path]

    #add other choosen parameters
    if api_key is not None:

        arguments += ['--api-key',api_key]

    if assembly_level is not None:
        arguments += ['--assembly-level',assembly_level]

    if reference and reference is not None:
        arguments += ['--reference']
    
    if exclude_atypical or exclude_atypical is None:
        arguments += ['--exclude-atypical']


    metadata = json.loads(subprocess.run(arguments,stdout=subprocess.PIPE,
                                         stderr=subprocess.PIPE,
                                         check=False).stdout)

    assemblies_ids = []

    if metadata['total_count'] != 0:
        #if any metadata with specified criteria were found
        metadata_all = metadata
        for meta in metadata_all['reports']:
            assemblies_ids.append(meta['accession'])

    else:
        #if no metadata was found, return json with {"total_count": 0}
        metadata_all = metadata

    assemblies_ids = []

    failed_list = [x for x in assemblies_ids if x not in input_ids]

    return [failed_list,metadata_all]

def metadata_from_id_list(id_list_path,size_threshold,max_contig_number,genome_size,
                          assembly_level,reference,verify_status,exclude_atypical,api_key,):
    """
    Function that from a list of ids and filtering criterea, filters the id list.

    Parameters
    ----------
    id_list_path:
    size_threshold: float
        (0 >= x >= 1).
    max_contig_number: int
        (>0).
    genome_size: int
        (>0).
    assembly_level: string
        containing assembly levels separated by "," .
    reference: Bool
    api_key: str

    Returns
    -------
    failed_list: list
        containing assemblies that failed criteria.
    accepted_list: list
        containing assemblies that passed criteria.
    """

    with open(id_list_path,'r',encoding='utf-8') as id_list:
        ids = list(csv.reader(id_list,delimiter='\t'))

    ids = [item for sublist in ids for item in sublist]

    if len(ids) == 0:
        os.sys.exit("\nInput file has no assemblies IDs.")

    #if input list has corrects format ids.
    if not all(i.startswith(('GCF_','GCA_')) for i in ids):
        sys.exit('\nOne or more ids in the input list have wrong format. Must '
                 'start with either GCF_ or GCA_')

    #if verification is needed
    if (genome_size is not None
        or size_threshold is not None
        or max_contig_number is not None
        or assembly_level is not None
        or reference and reference is not None):

        print("Verifying assemblies to be downloaded by specified criteria:")

        if genome_size is not None and size_threshold is not None:
            print(f"Genome size of: {genome_size}")
            print(f"Size threshold of: {size_threshold}")

        elif genome_size is None and size_threshold is None:
            print("Genome size of: Not specified")
            print("Size threshold of: Not specified")
        else:
            print("Both genome size and size threshold need to be specified.")
            print("Setting as:")
            print("    Genome size of: Not specified")
            print("    Size threshold of: Not specified")

        if max_contig_number is not None:
            print(f"Maximum number of contigs: {max_contig_number}")

        else:
            print("Maximum number of contigs: Not specified")

        if assembly_level is not None:
            print(f"assembly level at: {assembly_level}")

        else:
            print("assembly level at: Not specified (using Defaut: all)")

        if reference and reference is not None:
            print("Only reference genomes: True")

        else:
            print("Only reference genomes: False")

        if verify_status is not None:
            print(f"Remove suppressed assemblies: {verify_status}")
        else:
            print("Remove suppressed assemblies: True")
        
        if exclude_atypical is not None:
            print(f"Verify if assemblies are atypical: {exclude_atypical}")
        else:
            print("Verify if assemblies are atypical: True")
            
        verify_list = True

    else:
        verify_list = False
        print("No assembly verification criteria was specified...")
        print("Using default Criterias.")
        print("Proceeding...")

    accepted_list = []


    """
    Fetch metadata for all assemblies and verify according to desired fitering
    criteria
    """
    failed_list,metadata_all = metadata_fetcher_ids(ids,
                                                    id_list_path,
                                                    assembly_level,
                                                    reference,
                                                    exclude_atypical,
                                                    api_key)
    #Verify assemblies
    if metadata_all['total_count'] != 0:
        if verify_list:
            for metadata in metadata_all['reports']:
                if verify_assembly(metadata,
                                   size_threshold,
                                   max_contig_number,
                                   genome_size,
                                   verify_status):

                    accepted_list.append(metadata['accession'])
                else:
                    failed_list.append(metadata['accession'])

        else:
            accepted_list = ids

    return [failed_list,accepted_list]

def metadata_fetcher_taxon(taxon,assembly_level,reference,assembly_source,
                            exclude_atypical,api_key):
    """
    This function based on an input taxon fetches json object (dict).

    Parameters
    ----------
    taxon: string
    assembly_level: str
        containing assembly levels separated by "," .
    reference: Bool
    assembly_source: str
        'RefSeq' or 'GenBank'
    api_key: str

    Returns
    -------
    found_metadata: bool
    all_assemblies: list (containing all ids)
    metadata_filtered: json object (dict) for a all assemblies that
                        passed the criteria.
    """

    arguments = ['datasets','summary','genome','taxon', taxon]

    #filter by choosen assembly source
    if assembly_source is not None:
        arguments += ['--assembly-source',assembly_source]

    #find all possible assemblies ids
    metadata = json.loads(subprocess.run(arguments,stdout=subprocess.PIPE,
                                         stderr=subprocess.PIPE,
                                         check=False).stdout)

    try:
        metadata_all = metadata['reports']
        has_summary = True
    except KeyError:
        all_assemblies = []
        has_summary = False
    else:
        all_assemblies = []

        for meta in metadata_all:

            all_assemblies.append(meta['accession'])

    #add other choosen parameters
    if api_key is not None:

        arguments += ['--api-key',api_key]

    if assembly_level is not None:
        arguments += ['--assembly-level',assembly_level]

    if reference and reference is not None:
        arguments += ['--reference']
        
    if exclude_atypical or exclude_atypical is None:
        arguments += ['--exclude-atypical']
        
    if has_summary:
        # find all metadata that pass initial criterias
        metadata_filtered = json.loads(subprocess.run(arguments,stdout=subprocess.PIPE,
                                                      stderr=subprocess.PIPE,
                                                      check=False).stdout)
    else:
        #If assemblies of that a taxon don't have summary it may be cause due to
        #writing wrong taxon
        metadata_filtered = {'total_count': 0}

    return [all_assemblies,metadata_filtered]

def metadata_from_taxon(taxon,size_threshold,max_contig_number,genome_size,
                          assembly_level,reference,assembly_source,verify_status,
                          exclude_atypical,api_key):
    """
    Fetches the ids that pass the filtering criteria.

    Parameters
    ----------
    taxon: str
    size_threshold: float
        (0 >= x >= 1).
    max_contig_number: int
        (>0).
    genome_size: int
        (>0).
    assembly_level: string
        containing assembly levels separated by "," .
    reference: Bool
    assembly_source: string
        'RefSeq' or 'GenBank'
    api_key: str

    Returns
    -------
    failed_list: list
        containing assemblies that failed criteria
    accepted_list: list
        containing assemblies that passed criteria
    """

    if (assembly_source is not None
        or genome_size is not None
        or size_threshold is not None
        or max_contig_number is not None
        or assembly_level is not None
        or reference is not None
        or verify_status is not None):
        """
        if verification is needed
        """

        print("\nVerifying assemblies to be downloaded by specified criteria:")

        print(f"Taxon : {taxon}")

        if assembly_source is not None:

            if assembly_source != "all":
                print(f"Fetching assemblies from: {assembly_source}")

            else:
                print("Fetching assemblies from: RefSeq,GenBank")

        else:
            print("Fetching assemblies from: RefSeq,GenBank")

        if genome_size is not None and size_threshold is not None:
            print(f"Genome size of: {genome_size}")
            print(f"Size threshold of: {size_threshold}")

        elif genome_size is None and size_threshold is None:
            print("Genome size of: Not specified")
            print("Size threshold of: Not specified")
        else:
            print("Both genome size and size threshold need to be specified.")
            print("Setting as:")
            print("    Genome size of: Not specified")
            print("    Size threshold of: Not specified")

        if max_contig_number is not None:
            print(f"Maximum number of contigs: {max_contig_number}")

        else:
            print("Maximum number of contigs: Not specified")

        if assembly_level is not None:
            print(f"assembly level at: {assembly_level}")

        else:
            print("assembly level at: Not specified (using Defaut: all)")

        if reference and reference is not None:
            print("Only reference genomes: True")

        else:
            print("Only reference genomes: False")

        if verify_status is not None:
            print(f"Remove suppressed assemblies: {verify_status}")
        else:
            print("Remove suppressed assemblies: True")
            
        if exclude_atypical is not None:
            print(f"Verify if assemblies are atypical: {exclude_atypical}")
        else:
            print("Verify if assemblies are atypical: True")

        verify_list = True

    else:
        """
        If no verification is needed.
        """
        verify_list = False
        print("No assembly verification criteria was specified...")
        print("Using default Criterias.")
        print("Proceeding...")

    #get all possible ids and all ids filtered by assembly_level or reference
    all_ids,metadata_filtered = metadata_fetcher_taxon(taxon,
                                                        assembly_level,
                                                        reference,
                                                        assembly_source,
                                                        exclude_atypical,
                                                        api_key)
    accepted_list = []
    failed_list = []

    #verify by certain criteria
    if metadata_filtered['total_count'] != 0:
        if verify_list:
            for metadata in metadata_filtered['reports']:

                if verify_assembly(metadata,size_threshold,max_contig_number,
                                   genome_size,verify_status):

                    accepted_list.append(metadata['accession'])

                else:

                    failed_list.append(metadata['accession'])
        else:
            for metadata in metadata_filtered['reports']:

                accepted_list.append(metadata['accession'])

    #add ids that were initially removed by assembly_level or reference to failed
    #list
    failed_list = [x for x in all_ids if x not in accepted_list] + failed_list

    return [failed_list,accepted_list]
