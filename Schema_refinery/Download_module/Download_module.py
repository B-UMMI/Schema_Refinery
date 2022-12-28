#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec 15 18:18:38 2022

@author: mykyta
"""
import os
import csv
import subprocess
import ast

def tryeval(val):
  try:
    val = ast.literal_eval(val)
  except ValueError:
    pass
  return val

def main(args):
    
    print("Fetching assemblies from {}".format(args.database))    
    
    if args.filter_criteria_path != None:
        
        with open(args.filter_criteria_path,'r') as filters:
            criterias = dict(csv.reader(filters, delimiter='\t'))
        
        #Transform dictionary keys and values into variables
        abundance = tryeval(criterias['abundance'])
        genome_size = tryeval(criterias['genome_size'])
        size_threshold = tryeval(criterias['size_threshold'])
        max_contig_number = tryeval(criterias['max_contig_number'])
        known_st = tryeval(criterias['known_st'])
        any_quality = tryeval(criterias['any_quality'])
        ST_list_path = tryeval(criterias['ST_list_path'])
        assembly_level = tryeval(criterias['assembly_level'])
        reference = tryeval(criterias['reference_genome'])
        assembly_source = tryeval(criterias['assembly_source'])
        file_extension = tryeval(criterias['file_extension'])

           
    
    if args.database == 'NCBI':
        """
        If Database option is set to NCBI, run one of the two tools depending
        if there is a input table containing the ids of desired assemblies, if
        table exists then script downloads them otherwise it downloads all 
        assemblies for the desired species.
        """
        
        try:
            from Download_module import ncbi_datasets_summary
            from Download_module import ncbi_linked_ids
            from Download_module import fetch_metadata
                
        except ModuleNotFoundError:
            from Schema_refinery.Download_module import ncbi_datasets_summary
            from Schema_refinery.Download_module import ncbi_linked_ids
            from Schema_refinery.Download_module import fetch_metadata    
        
        if args.input_table is not None:
            
            failed_list,list_to_download = ncbi_datasets_summary.metadata_from_id_list(args.input_table,
                                                          size_threshold,
                                                          max_contig_number,
                                                          genome_size,
                                                          assembly_level,
                                                          reference,
                                                          args.threads,
                                                          args.api_key)
            
            #save ids to download
            with open(os.path.join(args.output_directory,
                                   "assemblies_ids_to_download.tsv"),'w+') as ids_to_txt:
                
                ids_to_txt.write("\n".join(map(str, list_to_download)))
            
            #save ids that failed criteria
            with open(os.path.join(args.output_directory,
                                   "id_failed_criteria.tsv"),'w+') as ids_to_txt:
                
                ids_to_txt.write("\n".join(map(str, failed_list)))                
            
            arguments = ['datasets','download','genome','accession','--inputfile',
                         os.path.join(args.output_directory,
                                      'assemblies_ids_to_download.tsv')]
            
            if args.api_key != None:
                arguments = arguments + ['--api-key',args.api_key]
            
            #download assemblies
            if args.download:
                os.chdir(args.output_directory)
                subprocess.run(arguments)

                
        else:
            
            failed_list,list_to_download = ncbi_datasets_summary.metadata_from_species(args.species,
                                                          size_threshold,
                                                          max_contig_number,
                                                          genome_size,
                                                          assembly_level,
                                                          reference,
                                                          args.api_key)
            
            
            #save ids to download
            with open(os.path.join(args.output_directory,
                                   "assemblies_ids_to_download.tsv"),'w+') as ids_to_txt:
                
                ids_to_txt.write("\n".join(map(str, list_to_download)))
            
            #save ids that failed criteria
            with open(os.path.join(args.output_directory,
                                   "id_failed_criteria.tsv"),'w+') as ids_to_txt: 
                
                ids_to_txt.write("\n".join(map(str, failed_list)))
 
            arguments = ['datasets','download','genome','taxon','--inputfile',
                         os.path.join(args.output_directory,
                                      'assemblies_ids_to_download.tsv')]
            
            if args.api_key != None:
                arguments = arguments + ['--api-key',args.api_key]
            
            if args.download:
                os.chdir(args.output_directory)
                subprocess.run(arguments)
        
        if args.f_metadata:
                
            if not os.path.exists(os.path.join(args.output_directory,'metadata')):
                os.mkdir(os.path.join(args.output_directory,'metadata'))  
                    
                
            print("\nFetching related ids...")
            ncbi_linked_ids.main(os.path.join(args.output_directory,'assemblies_ids_to_download.tsv'),
                                 os.path.join(args.output_directory,'metadata/id_matches.tsv'),
                                 args.email, 
                                 args.threads, 
                                 args.retry, 
                                 args.api_key)
                
            print("\nFetching additional metadata...")                
            fetch_metadata.main(os.path.join(args.output_directory,'assemblies_ids_to_download.tsv'),
                                os.path.join(args.output_directory,'metadata'),
                                args.email,
                                args.threads,
                                args.api_key,
                                args.retry)
            
    else:
        try:
            from Download_module import ena661k_assembly_fetcher
            from Download_module import ncbi_linked_ids
            from Download_module import fetch_metadata
            
        except ModuleNotFoundError:
            from Schema_refinery.Download_module import ena661k_assembly_fetcher
            from Schema_refinery.Download_module import ncbi_linked_ids
            from Schema_refinery.Download_module import fetch_metadata
            
        ena661k_assembly_fetcher.main(args.metadata_table,
                                      args.paths_table,
                                      args.species, 
                                      args.output_directory,
                                      args.download, 
                                      abundance, 
                                      genome_size, 
                                      size_threshold,
                                      max_contig_number, 
                                      args.species, 
                                      known_st, 
                                      any_quality, 
                                      args.stride,
                                      args.retry, 
                                      ST_list_path, 
                                      args.threads)

        if args.f_metadata:
            
            if not os.path.exists(os.path.join(args.output_directory,'metadata')):
                os.mkdir(os.path.join(args.output_directory,'metadata'))  
            
            ncbi_linked_ids.main(os.path.join(args.output_directory,'assemblies_ids.tsv'),
                                 os.path.join(args.output_directory,'metadata/id_matches.tsv'),
                                 args.email, 
                                 args.threads, 
                                 args.retry, 
                                 args.api_key)
            
            fetch_metadata.main(os.path.join(args.output_directory,'assemblies_ids.tsv'),
                                os.path.join(args.output_directory,'metadata'),
                                args.email,
                                args.threads,
                                args.api_key,
                                args.retry)
        