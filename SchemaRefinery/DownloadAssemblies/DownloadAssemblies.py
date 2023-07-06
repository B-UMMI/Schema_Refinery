#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
"""


import os
import sys
import subprocess
import pandas as pd

try:
    from DownloadAssemblies import ena661k_assembly_fetcher
    from DownloadAssemblies import ncbi_datasets_summary
    from DownloadAssemblies import ncbi_linked_ids
    from DownloadAssemblies import fetch_metadata
except ModuleNotFoundError:
    from SchemaRefinery.DownloadAssemblies import ena661k_assembly_fetcher
    from SchemaRefinery.DownloadAssemblies import ncbi_datasets_summary
    from SchemaRefinery.DownloadAssemblies import ncbi_linked_ids
    from SchemaRefinery.DownloadAssemblies import fetch_metadata

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
    conda_path = subprocess.Popen(['conda', 'info'], stdout=subprocess.PIPE)

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

    return os.path.join(sr_path, 'ena661k_files')

def get_all_metadata(metadata_dir, args):
    linked_ids_file = os.path.join(metadata_dir, 'id_matches.tsv')
    ids_file = os.path.join(metadata_dir, 'assemblies_ids_to_download.tsv')

    print("\nFetching RefSeq, Genbank and SRA IDs linked and BioSample ID...")
    ncbi_linked_ids.main(ids_file,
                            linked_ids_file,
                            args.email,
                            args.threads,
                            args.retry,
                            args.api_key)

    # fetch BioSample identifiers
    biosamples = pd.read_csv(linked_ids_file,
                                delimiter='\t')['BioSample'].values.tolist()
    # exclude samples without BioSample identifier
    biosamples = [i for i in biosamples if type(i) is str]
    # save BioSample identifiers to file
    biosample_file = os.path.join(metadata_dir, 'biosamples.tsv')
    with open(biosample_file, 'w+', encoding='utf-8') as ids:
        ids.write("\n".join(biosamples)+'\n')

    print("\nFetching metadata associated to the BioSample ID...")
    fetch_metadata.main(biosample_file,
                        metadata_dir,
                        args.email,
                        args.threads,
                        args.api_key,
                        args.retry)
    
    os.remove(biosample_file)

def main(args):
    # if both input table and taxon are present
    if args.input_table is not None and args.taxon is not None:
        sys.exit("\nError: Downloading from input table or downloading "
                 "by taxon are mutually exclusive.")

    # Table or taxon name must be present
    if args.input_table is None and args.taxon is None:
        sys.exit("\nError: Must provide an input table -i/--input-table or a taxon name -t/--taxon.")

    # if input table and database is set as ENA661k.
    if args.input_table is not None and 'ENA661K' in args.database:
        sys.exit("\nError: Only assemblies from NCBI can be fetched "
                 "from a input file. ENA661K was parsed.")

    if args.filtering_criteria:
        criteria = args.filtering_criteria
        
        if criteria['assembly_source'] == ['GenBank'] and (criteria['verify_status'] is True or criteria['verify_status'] is None):
            sys.exit("\nError: Assembly status can only be verified for assemblies obtained from RefSeq (Set to False, Default(None) = True)")
    else:
        criteria = None

    # create output directory
    if os.path.isdir(args.output_directory) is False:
        os.mkdir(args.output_directory)

    print(f"\nFetching assemblies from {args.database} datasets.")
    if 'NCBI' in args.database:
        # user provided NCBI Genome Assembly and Annotation report
        if args.input_table is not None:
            # read TSV file that contains NCBI report
            with open(args.input_table,'r') as id_list:
                assembly_ids = id_list.read().splitlines()

            if len(assembly_ids) == 0:
                print("No assembly identifiers provided.")

            if criteria is not None:
                metadata = ncbi_datasets_summary.fetch_metadata(args.input_table, None, criteria, args.api_key)

                if metadata['total_count'] == 0:
                    os.sys.exit("\nNo assemblies that satisfy the selected criteria were found.")
        else:
            # Fetch from taxon identifier
            metadata = ncbi_datasets_summary.fetch_metadata(None, args.taxon, criteria, args.api_key)
            
            if metadata['total_count'] == 0:
                os.sys.exit("\nNo assemblies that satisfy the selected criteria were found.")

            assembly_ids = [sample['accession'] for sample in metadata['reports']]

        total_ids = len(assembly_ids)

        failed = []
        passed = []
        if criteria is not None:
            # Validate assemblies
            if metadata['total_count'] > 0:
                for sample in metadata['reports']:
                    current_accession = sample['accession']
                    valid = ncbi_datasets_summary.verify_assembly(sample,
                                            criteria['size_threshold'],
                                            criteria['max_contig_number'],
                                            criteria['genome_size'],
                                            criteria['verify_status'])
                    if valid is True:
                        passed.append(current_accession)
                    else:
                        failed.append(current_accession)
                assembly_ids = passed

            print(f"\n{len(assembly_ids)} passed filtering criteria.")

        metadata_ncbi_directory = os.path.join(args.output_directory, 'metadata_ncbi')
        if not os.path.exists(metadata_ncbi_directory):
            os.mkdir(metadata_ncbi_directory)

        # save ids to download
        valid_ids_file = os.path.join(metadata_ncbi_directory,
                                      "assemblies_ids_to_download.tsv")
        with open(valid_ids_file, 'w', encoding='utf-8') as ids_to_txt:
            ids_to_txt.write("\n".join(assembly_ids)+'\n')

        # save ids that failed criteria
        failed_ids_file = os.path.join(metadata_ncbi_directory,
                                       "id_failed_criteria.tsv")
        with open(failed_ids_file, 'w', encoding='utf-8') as ids_to_txt:
            ids_to_txt.write("\n".join(failed)+'\n')

        # If any assembly passed filtering criteria
        if len(assembly_ids) == 0:
            print("\nNo assemblies meet the desired filtering criteria.")
            sys.exit("\nAssemblies that failed are in the following "
                     "TSV file: {}".format(failed_ids_file))
        # download assemblies
        if args.download:
            # Build initial arguments for the subprocess run of datasets
            arguments = ['datasets', 'download', 'genome',
                         'accession', '--inputfile', valid_ids_file]

            if args.api_key is not None:
                arguments.extend(['--api-key', args.api_key])

            if criteria['file_to_include'] is not None:
                if "genome" not in criteria['file_to_include']:
                    criteria['file_to_include'].append("genome")

                arguments.extend(['--include', ','.join(criteria['file_to_include'])])

            assemblies_zip = os.path.join(args.output_directory, 'assemblies.zip')
            arguments.extend(['--filename', assemblies_zip])
            print("\nDownloading assemblies...")
            subprocess.run(arguments, check=False)
        else:
            print("\nThe list of identifiers for the assemblies that passed the filtering criteria was saved to: {}".format(valid_ids_file))

    # Download from ENA661k
    if 'ENA661K' in args.database:
        # Path for ena661k files
        try:
            # Find local conda path
            sr_path = find_local_conda_env()
        except:
            # if not using conda, use output directory instead
            sr_path = os.path.join(args.output_directory, 'ena661k_files')

        metadata_ena_directory = ena661k_assembly_fetcher.main(sr_path,
                                                           args.taxon,
                                                           args.output_directory,
                                                           args.download,
                                                           criteria,
                                                           args.retry,
                                                           args.threads)

    if args.fetch_metadata:
        if 'NCBI' in args.database:
            print("\nFetching metadata for NCBI assemblies...")
            get_all_metadata(metadata_ncbi_directory, args)
        if 'ENA661K' in args.database:
            print("\nFetching metadata for ENA661K assemblies...")
            get_all_metadata(metadata_ena_directory, args)
