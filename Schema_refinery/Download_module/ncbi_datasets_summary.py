#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Dec 27 13:05:38 2022

@author: mykyta
"""

import subprocess
import json
import concurrent.futures
import csv
from itertools import repeat

def verify_assembly(metadata_assembly,size_threshold,max_contig_number,
                    genome_size):
    """
    This function verifies assemblies by certain inputa criteria.
        
        Input: 
            metadata_assembly: json object (dict) for a single assembly
            size_threshold: float (0 >= x >= 1)
            max_contig_number: int (>0)
            genome_size: int (>0)
            
        Output:
            return : Boolean value (in order to see if passed or failed)
    """
    
    if genome_size is not None and size_threshold is not None:
        
        bot_limit = genome_size - (genome_size*size_threshold)
        top_limit = genome_size + (genome_size*size_threshold)
        
        if int(metadata_assembly['total_sequence_length']) >= top_limit:
            
            return False
        
        if int(metadata_assembly['total_sequence_length']) <= bot_limit:
            
            return False
    
    if max_contig_number is not None:
        
        if metadata_assembly['number_of_contigs'] >= max_contig_number:
            
            return False

    return True

                

def metadata_fetcher_id(input_id,assembly_level,reference,api_key):
    """
    This function based on an input id fetches json object (dict).
    
        Input:
            input_id: string (starts with GCF_ or GCA_)
            assembly_level: string (containing assembly levels separated by ",")
            reference: Boolean value
            api_key: string
            
        Output:
            return: 
                Input_id: string (starts with GCF_ or GCA_)
                metadata_assembly: json object (dict) for a single assembly
    """
    
    arguments = ['datasets','summary','genome','accession',input_id]
    
    if api_key is not None:
        
        arguments += ['--api-key',api_key]
    
    if assembly_level is not None:
        arguments += ['--assembly-level',assembly_level]
    
    if reference and reference is not None:
        arguments += ['--reference']
        
    try:
        metadata = subprocess.run(arguments,stdout=subprocess.PIPE)
                                   
        metadata_assembly = json.loads(metadata.stdout)['reports'][0]['assembly_stats']         
                  
    except:
        return [input_id,'Failed']

    return [input_id,metadata_assembly]

def metadata_from_id_list(id_list_path,size_threshold,max_contig_number,genome_size,
                          assembly_level,reference,threads,api_key):
    """
    Function that from a list of ids and filtering criterea, filters the id list.
    
        Input:
            id_list_path:
            size_threshold: float (0 >= x >= 1)
            max_contig_number: int (>0)
            genome_size: int (>0)
            assembly_level: string (containing assembly levels separated by ",")
            reference: Boolean value
            threads: int (>0)
            api_key: string
        
        Output:
            return:
                failed_list: list (containing assemblies that failed criteria)
                accepted_list: list (containing assemblies that passed criteria)
    """

    with open(id_list_path,'r') as id_list:
        ids = list(csv.reader(id_list,delimiter='\t'))

    ids = [item for sublist in ids for item in sublist]
    
    #if verification is needed
    if (genome_size is not None
        or size_threshold is not None 
        or max_contig_number is not None 
        or assembly_level is not None
        or reference and reference is not None):
        
        print("Verifying assemblies to be downloaded by specified criteria:")
        
        if genome_size is not None and size_threshold is not None:
            print("Genome size of: {}".format(genome_size))
            print("Size threshold of: {}".format(size_threshold))
            
        else:
            print("Genome size of: Not specified")
            print("Size threshold of: Not specified")
            
        if max_contig_number is not None:
            print("Maximum number of contigs: {}".format(max_contig_number))
            
        else:
            print("Maximum number of contigs: Not specified")
            
        if assembly_level is not None:
            print("assembly level at: {}".format(assembly_level))
            
        else:
            print("assembly level at: Not specified (using Defaut: all)")
            
        if reference and reference is not None:
            print("Only reference genomes: True")
            
        else:
            print("Only reference genomes: False")
        
        verify_list = True
    
    else:
        verify_list = False
        print("No assembly verification criteria was specified...")
        print("Using default Criterias.")
        print("Proceeding...")
    
    accepted_list = []
    failed_list = []
    
    #Verify list of ids
    if verify_list:
        with concurrent.futures.ThreadPoolExecutor(max_workers=threads) as executor:
    
            for metadata_assembly in executor.map(metadata_fetcher_id, 
                                                  ids,
                                                  repeat(assembly_level),
                                                  repeat(reference),
                                                  repeat(api_key)):
                
                if 'Failed' not in metadata_assembly:
                    
                    if verify_assembly(metadata_assembly[1],
                                           size_threshold,
                                           max_contig_number,genome_size):
                        
                        accepted_list.append(metadata_assembly[0])
                else:
                    failed_list.append(metadata_assembly[0])
        
    else:
        accepted_list = ids

    return [failed_list,accepted_list]

def metadata_fetcher_specie(species,assembly_level,reference,assembly_source,
                            api_key):
    """
    This function based on an input species fetches json object (dict).
    
        Input:
            species: string
            assembly_level: string (containing assembly levels separated by ",")
            reference: Boolean value
            assembly_source: string (RefSeq or GenBank)
            api_key: string
            
        Output:
            all_assemblies: list (containing all ids)
            metadata_filtered: json object (dict) for a all assemblies that
                                passed the criteria.
    """
    
    arguments = ['datasets','summary','genome','taxon', species]
    
    #find all possible assemblies ids
    metadata = subprocess.run(arguments,stdout=subprocess.PIPE)
                                   
    metadata_all = json.loads(metadata.stdout)['reports']
    
    all_assemblies = []
    
    for m in metadata_all:
        
        all_assemblies.append(m['accession'])
    
    if api_key is not None:
        
        arguments += ['--api-key',api_key]
    
    if assembly_level is not None:
        arguments += ['--assembly-level',assembly_level]
    
    if reference and reference is not None:
        arguments += ['--reference']
        
    if assembly_source is not None:
        arguments += ['--assembly-source',assembly_source]
        
    # find all metadata that pass initial criterias
    metadata = subprocess.run(arguments,stdout=subprocess.PIPE)
                                   
    metadata_filtered = json.loads(metadata.stdout)['reports']      
    
    return [all_assemblies,metadata_filtered]
    
def metadata_from_species(species,size_threshold,max_contig_number,genome_size,
                          assembly_level,reference,assembly_source,api_key):
    """
    Fetches the ids that pass the filtering criteria.
    
        Input:
            species: string
            size_threshold: float (0 >= x >= 1)
            max_contig_number: int (>0)
            genome_size: int (>0)
            assembly_level: string (containing assembly levels separated by ",")
            reference: Boolean value
            assembly_source: string (RefSeq or GenBank)
            api_key: string
            
        Output:
            return:
                failed_list: list (containing assemblies that failed criteria)
                accepted_list: list (containing assemblies that passed criteria)
    """
    #if verification is needed
    if (assembly_source is not None
        or genome_size is not None
        or size_threshold is not None 
        or max_contig_number is not None
        or assembly_level is not None
        or reference is not None):
        
        print("\nVerifying assemblies to be downloaded by specified criteria:")
        
        if assembly_source is not None:
            
            if assembly_source != "all":
                print("Fetching assemblies from: {}".format(assembly_source))
            
            else:
                print("Fetching assemblies from: RefSeq,GenBank")
            
        else:
            print("Fetching assemblies from: RefSeq,GenBank")
        
        if genome_size is not None and size_threshold is not None:
            print("Genome size of: {}".format(genome_size))
            print("Size threshold of: {}".format(size_threshold))
            
        else:
            print("Genome size of: Not specified")
            print("Size threshold of: Not specified")
            
        if max_contig_number is not None:
            print("Maximum number of contigs: {}".format(max_contig_number))
            
        else:
            print("Maximum number of contigs: Not specified")
            
        if assembly_level is not None:
            print("assembly level at: {}".format(assembly_level))
            
        else:
            print("assembly level at: Not specified (using Defaut: all)")
            
        if reference and reference is not None:
            print("Only reference genomes: True")
            
        else:
            print("Only reference genomes: False")
        
        verify_list = True
    
    else:
        verify_list = False
        print("No assembly verification criteria was specified...")
        print("Using default Criterias.")
        print("Proceeding...")

    #get all possible ids and all ids filtered by assembly_level or reference
    all_ids,metadata_filtered = metadata_fetcher_specie(species,
                                                        assembly_level,
                                                        reference,
                                                        assembly_source,
                                                        api_key)
    
    accepted_list = []
    failed_list = []
    
    #verify by certain criteria
    if verify_list:
        
        for metadata in metadata_filtered:
            
            if verify_assembly(metadata['assembly_stats'],size_threshold,
                               max_contig_number,genome_size):
                
                accepted_list.append(metadata['accession'])
            
            else:
                
                failed_list.append(metadata['accession'])
    else:
        
        for metadata in metadata_filtered:
            
            accepted_list.append(metadata['accession'])
    
    #add ids that were initially removed by assembly_level or reference to failed
    #list
    failed_list = [x for x in all_ids if x not in accepted_list] + failed_list
            
    return [failed_list,accepted_list]
    
                
                