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
    
    if genome_size is not None and size_threshold is not None:
        
        bot_limit = genome_size - (genome_size*size_threshold)
        top_limit = genome_size + (genome_size*size_threshold)
        
        if metadata_assembly['total_sequence_length'] > top_limit:
            
            return False
        
        if metadata_assembly['total_sequence_length'] < bot_limit:
            
            return False
    
    if max_contig_number is not None:
        
        if metadata_assembly['number_of_contigs'] > max_contig_number:
            
            return False
        
    else:
        
        return True

                

def metadata_fetcher_id(input_id,assembly_level,reference,api_key):
    
    arguments = ['datasets','summary','genome','accession',input_id]
    
    if api_key != None:
        
        arguments = arguments + ['--api-key',api_key]
    
    if assembly_level != None:
        arguments = arguments + ['--assembly-level',assembly_level]
    
    if reference and reference != None:
        arguments = arguments + ['--reference']
        
    try:
        metadata = subprocess.run(arguments,stdout=subprocess.PIPE)
                                   
        metadata_assembly = json.loads(metadata.stdout)['reports'][0]['assembly_stats']         
                  
    except:
        return [input_id,'Failed']

    return [input_id,metadata_assembly]

def metadata_from_id_list(id_list_path,size_threshold,max_contig_number,genome_size,
                          assembly_level,reference,threads,api_key):

    with open(id_list_path,'r') as id_list:
        ids = list(csv.reader(id_list,delimiter='\t'))

    ids = [item for sublist in ids for item in sublist]
    
    #if verification is needed
    if (genome_size is not None
        or size_threshold is not None 
        or max_contig_number is not None 
        or assembly_level != None
        or reference and reference != None):
        
        print("Verifying assemblies to be downloaded by specified criteria")
        
        if genome_size is not None and size_threshold is not None:
            print("\nGenome size of: {}".format(genome_size))
            print("\and size threshold of: {}".format(size_threshold))
            
        if max_contig_number is not None:
            print("Maximum number of contigs: {}".format(max_contig_number))
            
        if assembly_level != None:
            print("assembly level at: {}".format(assembly_level))
            
        if reference and reference != None:
            print("Only reference genomes: True")
        
        verify_list = True
    
    else:
        
        verify_list = False
    
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

def metadata_fetcher_specie(species,assembly_level,reference,api_key):
    
    arguments = ['datasets','summary','genome','taxon',species]
    
    #find all possible assemblies ids
    metadata = subprocess.run(arguments,stdout=subprocess.PIPE)
                                   
    metadata_all = json.loads(metadata.stdout)['reports']
    
    all_assemblies = []
    
    for m in metadata_all:
        
        all_assemblies.append(m['accession'])
    
    if api_key != None:
        
        arguments = arguments + ['--api-key',api_key]
    
    if assembly_level != None:
        arguments = arguments + ['--assembly-level',assembly_level]
    
    if reference and reference != None:
        arguments = arguments + ['--reference']
        

    metadata = subprocess.run(arguments,stdout=subprocess.PIPE)
                                   
    metadata_filtered = json.loads(metadata.stdout)['reports']      
    
    return [all_assemblies,metadata_filtered]
    
def metadata_from_species(species,size_threshold,max_contig_number,genome_size,
                          assembly_level,reference,api_key):
    
    #if verification is needed
    if (genome_size is not None
        or size_threshold is not None 
        or max_contig_number is not None):
        
        print("Verifying assemblies to be downloaded by specified criteria")
        
        if genome_size is not None and size_threshold is not None:
            print("\nGenome size of: {}".format(genome_size))
            print("\and size threshold of: {}".format(size_threshold))
            
        if max_contig_number is not None:
            print("Maximum number of contigs: {}".format(max_contig_number))
            
        if assembly_level != None:
            print("assembly level at: {}".format(assembly_level))
            
        if reference and reference != None:
            print("Only reference genomes: True")
        
        verify_list = True
    
    else:
        
        verify_list = False

    #get all possible ids and all ids filtered by assembly_level or reference
    all_ids,metadata_filtered = metadata_fetcher_specie(species,
                                                        assembly_level,
                                                        reference,
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
    
    #add ids that were initially removed
    failed_list = [x for x in all_ids if x not in accepted_list] + failed_list
            
    return [failed_list,accepted_list]
    
                
                