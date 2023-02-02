#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec 15 18:18:38 2022

AUTHOR

    Mykyta Forofontov
    github: @MForofontov
"""
import os
import sys
import csv
import subprocess
import ast
import pandas as pd

try:
    from DownloadAssemblies import ena661k_assembly_fetcher
    from DownloadAssemblies import ncbi_datasets_summary
    from DownloadAssemblies import ncbi_linked_ids
    from DownloadAssemblies import fetch_metadata

except ModuleNotFoundError:
    from Schema_refinery.DownloadAssemblies import ena661k_assembly_fetcher
    from Schema_refinery.DownloadAssemblies import ncbi_datasets_summary
    from Schema_refinery.DownloadAssemblies import ncbi_linked_ids
    from Schema_refinery.DownloadAssemblies import fetch_metadata

def tryeval(val):
    """
    Evaluates the type of the input.

    Parameter
    ---------
    val : any

    Returns
    -------
    val : any
        converted to the right type.
    """

    try:
        val = ast.literal_eval(val)
    except ValueError:
        pass
    return val

def find_local_conda_env():
    """
    This function fetches current conda env path by using 'conda info'.

    Parameter
    ---------
    None

    Returns
    -------
    return: str
        containing path to folder inside schema refinery
        env, may not exist if this module was not used to download
        assemblies from ENA661K.´
    """
    conda_path = subprocess.Popen(['conda','info'],stdout=subprocess.PIPE)

    stdout, stderr = conda_path.communicate()

    conda_info = stdout.decode("utf-8").splitlines()

    info_dict = {}

    for info in conda_info:

        if info == '':
            continue

        key_value = info.strip().split(':')

        if len(key_value) == 2:
            info_dict[key_value[0].strip()] = [key_value[1].strip()]

            out_key = key_value[0].strip()
        else:
            info_dict[out_key].append(key_value[0].strip())

    sr_path = info_dict['active env location'][0]

    return os.path.join(sr_path,'ena661k_files')

def filtering_criteria_variables(filtering_criteria_path, expected_criterias):
    """
    This function imports and verifies filtering criteria obtained from an input
    table.

    Parameter
    ---------
    filtering_criteria_path : str
        filtering criteria path.
    expected_criterias : list
        with expected criterias.

    Returns
    -------
    return: 11 filtering criteria variables, with their type and format
            verified, some may be None values.
    """

    with open(filtering_criteria_path,'r',encoding='utf-8') as filters:
        criterias = dict(csv.reader(filters, delimiter='\t'))

    for key, value in criterias.items():
        criterias[key] = tryeval(value)

    unexpected_keys = [x for x in criterias if x not in expected_criterias]

    exit_bool = False

    if len(unexpected_keys) > 0:
        print("\nError: Following unexpected parameters:")
        for k in unexpected_keys:
            print(k)

        exit_bool = True

    missing_keys = [x for x in expected_criterias if x not in criterias]

    if len(missing_keys) > 0:
        print("\nError: Missing following parameters:")
        for k in missing_keys:
            print(k)

        exit_bool = True

    if exit_bool:
        os.sys.exit()


    print("\nLoaded filtering criteria file successfully.")

    #Verify variables integrity and validity
    print("\nVerifying filtering criteria input table.")

    wrong_inputs = []

    if criterias['abundance'] is not None:
        if not 0 < criterias['abundance'] <= 1:
            wrong_inputs.append('abundance: float values between 0 < x <= 1')

    if criterias['genome_size'] is not None:
        if (criterias['genome_size'] < 0 or type(criterias['genome_size']) is not int):
            wrong_inputs.append('genome_size: int value > 0')

    if criterias['size_threshold'] is not None:
        if not 0 < criterias['size_threshold'] <= 1:
            wrong_inputs.append('size_threshold: float values between 0 < x <= 1')

    if criterias['max_contig_number'] is not None:
        if (criterias['max_contig_number'] < 0 or type(criterias['max_contig_number']) is not int):
            wrong_inputs.append('max_contig_number: int value > 0')

    if criterias['known_st']  is not None:
        if type(criterias['known_st']) is not bool:
            wrong_inputs.append('known_st: bool value')

    if criterias['any_quality'] is not None:
        if type(criterias['any_quality']) is not bool:
            wrong_inputs.append('any_quality: bool value')

    if criterias['ST_list_path'] is not None:
        if (not os.path.exists(criterias['ST_list_path'])
            or type(criterias['ST_list_path']) is not str):
            wrong_inputs.append('ST_list_path: path to the file with ST list')

    if criterias['assembly_level'] is not None:
        if type(criterias['assembly_level']) is str:
            if (not all(level in ['chromosome','complete','contig','scaffold']
                    for level in criterias['assembly_level'].split(','))):

                wrong_inputs.append('assembly_level: '
                                    'one or more of the following separated '
                                    'by a comma: '
                                    'chromosome,complete,contig,scaffold')
        else:
            wrong_inputs.append('assembly_level: '
                                'one ore more of the following separated '
                                'by a comma: '
                                'chromosome,complete,contig,scaffold')

    if criterias['reference_genome'] is not None:
        if type(criterias['reference_genome']) is not bool:
            wrong_inputs.append('reference: bool value')

    if criterias['assembly_source'] is not None:
        if (criterias['assembly_source'] not in ['RefSeq','GenBank','all']):
            wrong_inputs.append('assembly_source: one of the following:'
                                'RefSeq,GenBank,all')

    if criterias['file_extension']  is not None:
        if (criterias['file_extension']  not in ['genome','rna','protein','cds','gff3','gtf',
                                                 'gbff','seq-report','none']):
            wrong_inputs.append('file_extension: one or more of the following separated '
                                'by comma: genome,rna,protein,cds,gff3,gtf, '
                                'gbff,seq-report,none')

    if criterias['verify_status'] is not None:
        if type(criterias['verify_status']) is not bool:
            wrong_inputs.append('verify_status: bool or None')

    if len(wrong_inputs) > 0:
        print('\nFiltering table inputs have wrong values or types:')
        for inputs in wrong_inputs:
            print('\n' + inputs)
        os.sys.exit()

    else:
        return criterias

def main(args):
    """
    Main function that is called at schema_refinery.

    From the input arguments this function determines which function to use to
    download assemblies.
    """
    if args.input_table is not None and args.database == 'ENA661K':
        sys.exit("\nError: Only assemblies from NCBI can be fetched from a input file.")


    expected_criterias = ['abundance', 'genome_size', 'size_threshold', 'max_contig_number',
                          'known_st', 'any_quality', 'ST_list_path', 'assembly_level',
                          'reference_genome', 'assembly_source', 'file_extension',
                          'verify_status']

    criterias = {}

    if args.filter_criteria_path is not None:
        #import filtering criterias file
        criterias = filtering_criteria_variables(args.filter_criteria_path,
                                                 expected_criterias)

    else:
        for criteria in expected_criterias:
            criterias[criteria] = None

    print("\nFetching assemblies from {} datasets.".format(args.database))

    if 'NCBI' in args.database:
        """
        If Database option is set to NCBI, run one of the two tools depending
        if there is a input table containing the ids of desired assemblies, if
        table exists then script downloads them otherwise it downloads all
        assemblies for the desired species.
        """
        if args.input_table is not None:
            """
            If input table with ids is present.
            Based on the table, executes ncbi_datasets_summary.metadata_from_id_list
            function to obtain ids of assemblies that fail and pass the
            filtering criteria creating 2 output TSV files in the output folder.

            Based on ids that passed the filtering criterias, if --download is
            triggered download will be executed by running a subprocess of
            datasets.
            """
            failed_list,list_to_download = ncbi_datasets_summary.metadata_from_id_list(args.input_table,
                                                          criterias['size_threshold'],
                                                          criterias['max_contig_number'],
                                                          criterias['genome_size'],
                                                          criterias['assembly_level'],
                                                          criterias['reference_genome'],
                                                          criterias['verify_status'],
                                                          args.api_key)

            if not os.path.exists(os.path.join(args.output_directory,'metadata_ncbi')):
                os.mkdir(os.path.join(args.output_directory,'metadata_ncbi'))
            #save ids to download
            with open(os.path.join(args.output_directory,
                                   "metadata_ncbi/assemblies_ids_to_download.tsv"),
                      'w+', encoding='utf-8') as ids_to_txt:

                ids_to_txt.write("\n".join(map(str, list_to_download)))

            #save ids that failed criteria
            with open(os.path.join(args.output_directory,
                                   "metadata_ncbi/id_failed_criteria.tsv"),
                      'w+',encoding='utf-8') as ids_to_txt:

                ids_to_txt.write("\n".join(map(str, failed_list)))


            #If any assembly passed filtering criteria
            if not list_to_download:
                print("\nNo assemblies meet the desired filtering criterias.")

                sys.exit("\nAssemblies that failed are in the following"
                         "TSV file: {}".format(os.path.join(args.output_directory,
                                                            "metadata_ncbi/id_failed_criteria.tsv.tsv")))

            #Build initial arguments for the subprocess run of datasets
            arguments = ['datasets','download','genome','accession','--inputfile',
                         os.path.join(args.output_directory,
                                      'metadata_ncbi/assemblies_ids_to_download.tsv')]

            if args.api_key is not None:
                arguments += ['--api-key',args.api_key]

            if criterias['file_extension'] is not None:
                arguments += ['--include',criterias['file_extension']]

            #download assemblies
            if args.download:
                print("\nDownloading assemblies...")
                os.chdir(args.output_directory)
                subprocess.run(arguments, check=False)

            else:
                print("\nAssemblies to be downloaded are in the following TSV file: {}".format(
                    os.path.join(args.output_directory,"metadata_ncbi/assemblies_ids_to_download.tsv")))

        else:
            """
            If no table is present.
            Executes ncbi_datasets_summary.metadata_from_species function to
            obtain two list of assemblies that passed or failed the filtering
            criteria creating 2 output TSV files in the output folder.

            Based on ids that passed the filtering criterias, if --download is
            triggered download will be executed by running a subprocess of
            datasets.
            """
            #Fetch from species identifier
            failed_list,list_to_download = ncbi_datasets_summary.metadata_from_species(args.species,
                                                          criterias['size_threshold'],
                                                          criterias['max_contig_number'],
                                                          criterias['genome_size'],
                                                          criterias['assembly_level'],
                                                          criterias['reference_genome'],
                                                          criterias['assembly_source'],
                                                          criterias['verify_status'],
                                                          args.api_key)

            if not os.path.exists(os.path.join(args.output_directory,'metadata_ncbi')):
                os.mkdir(os.path.join(args.output_directory,'metadata_ncbi'))

            #save ids to download
            with open(os.path.join(args.output_directory,
                                   "metadata_ncbi/assemblies_ids_to_download.tsv"),'w+') as ids_to_txt:

                ids_to_txt.write("\n".join(map(str, list_to_download)))

            #save ids that failed criteria
            with open(os.path.join(args.output_directory,
                                   "metadata_ncbi/id_failed_criteria.tsv"),'w+') as ids_to_txt:

                ids_to_txt.write("\n".join(map(str, failed_list)))

            #If any assembly passed filtering criteria
            if not list_to_download:
                print("\nNo assemblies meet the desired filtering criterias.")

                sys.exit("\nAssemblies that failed are in the following "
                         "TSV file: {}".format(os.path.join(args.output_directory,
                                                            "metadata_ncbi/id_failed_criteria.tsv.tsv")))

            #Build initial arguments for the subprocess run of datasets
            arguments = ['datasets','download','genome','accession','--inputfile',
                         os.path.join(args.output_directory,
                                      'metadata_ncbi/assemblies_ids_to_download.tsv')]

            if args.api_key is not None:
                arguments = arguments + ['--api-key',args.api_key]

            if args.download:
                print("\nDownloading assemblies...")
                os.chdir(args.output_directory)
                subprocess.run(arguments, check=False)

            else:
                print("\nAssemblies to be downloaded are in the following TSV file: {}".format(
                    os.path.join(args.output_directory,"assemblies_ids_to_download.tsv")))

        if args.f_metadata:
            """
            If -fm is not triggered, this block of code executes two functions,
            ncbi_linked_ids.main and fetch_metadata.main, the first one associates
            the input ids with other SRA, Biosample ids, while the second one
            based on Biosample id obtainded previously fetches metadata related to
            that Biosample.
            """


            print("\nFetching related ids...")

            ncbi_linked_ids.main(os.path.join(args.output_directory,'metadata_ncbi/assemblies_ids_to_download.tsv'),
                                 os.path.join(args.output_directory,'metadata_ncbi/id_matches.tsv'),
                                 args.email,
                                 args.threads,
                                 args.retry,
                                 args.api_key)

            print("\nFetching additional metadata...")

            biosamples = pd.read_csv(os.path.join(args.output_directory,'metadata_ncbi/id_matches.tsv'),
                        delimiter = '\t')['BioSample'].values.tolist()

            #create file with biosamples ids
            with open(os.path.join(args.output_directory,
                                   'metadata_ncbi/biosamples.tsv'),
                      'w+', encoding='utf-8') as ids:

                ids.write("\n".join(map(str, biosamples)))

            fetch_metadata.main(os.path.join(args.output_directory,'metadata_ncbi/biosamples.tsv'),
                                os.path.join(args.output_directory,'metadata_ncbi'),
                                args.email,
                                args.threads,
                                args.api_key,
                                args.retry)

            #remove biosamples file
            os.remove(os.path.join(args.output_directory,'metadata_ncbi/biosamples.tsv'))
    if 'ENA661K' in args.database:
        """
        This block of code executes everything related to ENA661K assemblies
        download.
        Apart from downloading assemblies if --download  is triggered, it also
        creates two TSV files at output directory, one containing ids of assemblies
        that failes filtering criteria other file that passed.
        """
        #Path for ena661k files
        try:
            #Find local conda path
            sr_path = find_local_conda_env()
        except:
            #if not using conda, use output directory instead
            sr_path = os.path.join(args.output_directory, 'ena661k_files')

        ena661k_assembly_fetcher.main(sr_path,
                                      args.species,
                                      args.output_directory,
                                      args.download,
                                      criterias['abundance'],
                                      criterias['genome_size'],
                                      criterias['size_threshold'],
                                      criterias['max_contig_number'],
                                      criterias['known_st'],
                                      criterias['any_quality'],
                                      args.stride,
                                      args.retry,
                                      criterias['ST_list_path'],
                                      args.threads)

        if args.f_metadata:
            """
            If -fm is not triggered, this block of code executes two functions,
            ncbi_linked_ids.main and fetch_metadata.main, the first one associates
            the input ids with other SRA, Biosample ids, while the second one
            based on Biosample id fetches metadata related to that Biosample.
            """
            if not os.path.exists(os.path.join(args.output_directory,'metadata_ena661k')):
                os.mkdir(os.path.join(args.output_directory,'metadata_ena661k'))

            print("\nFetching related ids...")

            ncbi_linked_ids.main(os.path.join(args.output_directory,'metadata_ena661k/assemblies_ids_to_download.tsv'),
                                 os.path.join(args.output_directory,'metadata_ena661k/id_matches.tsv'),
                                 args.email,
                                 args.threads,
                                 args.retry,
                                 args.api_key)

            print("\nFetching additional metadata...")

            fetch_metadata.main(os.path.join(args.output_directory,'metadata_ena661k/assemblies_ids_to_download.tsv'),
                                os.path.join(args.output_directory,'metadata_ena661k'),
                                args.email,
                                args.threads,
                                args.api_key,
                                args.retry)
        