#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec 15 18:18:38 2022

@author: mykyta
"""
import os

def main(args):
    
    email = args.email
    metadata = args.f_metadata
    api_key = args.api_key
    
    print("Fetching assemblies from {}\n".format(args.database))
    
    if args.database == 'NCBI':
        """
        If Database option is set to NCBI, run one of the two tools depending
        if there is a input table containing the ids of desired assemblies, if
        table exists then script downloads them otherwise it downloads all 
        assemblies for the desired species.
        """
    
        if args.input_table is not None:
        
            try:
                from Download_module import ncbi_assembly_fetcher
                from Download_module import ncbi_linked_ids
                from Download_module import fetch_metadata
                
            except ModuleNotFoundError:
                from Schema_refinery.Download_module import ncbi_assembly_fetcher
                from Schema_refinery.Download_module import ncbi_linked_ids
                from Schema_refinery.Download_module import fetch_metadata
            
            ncbi_assembly_fetcher.main(args.input_table, 
                                       args.output_directory,
                                       args.file_extension, 
                                       args.ftp, 
                                       args.threads, 
                                       args.species, 
                                       args.retry)
        
            if metadata:
                
                ncbi_linked_ids.main(os.path.join(args.output_directory,'assemblies_ids_ncbi.tsv'),
                                     os.path.join(args.output_directory,'metadata/id_matches.tsv'),
                                     email, 
                                     args.threads, 
                                     args.retry, 
                                     api_key)
                
                fetch_metadata.main(os.path.join(args.output_directory,'assemblies_ids_ncbi.tsv'),
                                    os.path.join(args.output_directory,'metadata'),
                                    email,
                                    args.threads,
                                    api_key,
                                    args.retry)
                
        else:
            
            try:
                from Download_module import ncbi_assembly_fetcher
                from Download_module import ncbi_linked_ids
                from Download_module import fetch_metadata
                
            except ModuleNotFoundError:
                from Schema_refinery.Download_module import ncbi_assembly_fetcher
                from Schema_refinery.Download_module import ncbi_linked_ids
                from Schema_refinery.Download_module import fetch_metadata
                
            pass
            
    else:
            
        del args.database
        del args.email
        del args.f_metadata
        del args.api_key
        del args.input_table
        
        print("omitting optional arguments for NCBI Database")
        
        del args.file_extension
        del args.ftp
            
        try:
            from Download_module import ena661k_assembly_fetcher
            from Download_module import ncbi_linked_ids
            from Download_module import fetch_metadata
            
        except ModuleNotFoundError:
            from Schema_refinery.Download_module import ena661k_assembly_fetcher
            from Schema_refinery.Download_module import ncbi_linked_ids
            from Schema_refinery.Download_module import fetch_metadata
            
        ena661k_assembly_fetcher.main(**vars(args))

        if metadata:
            
            if not os.path.exists(os.path.join(args.output_directory,'metadata')):
                os.mkdir(os.path.join(args.output_directory,'metadata'))  
            
            ncbi_linked_ids.main(os.path.join(args.output_directory,'assemblies_ids.tsv'),
                                 os.path.join(args.output_directory,'metadata/id_matches.tsv'),
                                 email, 
                                 args.threads, 
                                 args.retry, 
                                 api_key)
            
            fetch_metadata.main(os.path.join(args.output_directory,'assemblies_ids.tsv'),
                                os.path.join(args.output_directory,'metadata'),
                                email,
                                args.threads,
                                api_key,
                                args.retry)
        