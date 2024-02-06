#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Purpose
-------
This file call all other modules of the packages, containing several arguments
as input depending on the desired module.

Code documentation
------------------
"""

import sys
import argparse

try:
    from DownloadAssemblies import DownloadAssemblies
    from utils import parameter_validation as pv
    from RefineSchema import SpuriousLociIdentification
except ModuleNotFoundError:
    from SchemaRefinery.DownloadAssemblies import DownloadAssemblies
    from SchemaRefinery.utils import parameter_validation as pv
    from SchemaRefinery.RefineSchema import SpuriousLociIdentification


def download_assemblies():

    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument('DownloadAssemblies', nargs='+',
                        help='')

    # Common arguments between databases
    parser.add_argument('-db', '--database', type=str,
                        required=True, dest='database',
                        nargs='+',
                        choices=['NCBI', 'ENA661K'],
                        help='Databases from which assemblies will '
                             'be downloaded.')

    parser.add_argument('-o', '--output-directory', type=str,
                        required=True, dest='output_directory',
                        help='Path to the output directory.')

    parser.add_argument('-e', '--email', type=str,
                        required=True, dest='email',
                        help='Email provided to Entrez.')

    parser.add_argument('-t', '--taxon', type=str,
                        required=False, dest='taxon',
                        help='Scientific name of the taxon.')

    parser.add_argument('-th', '--threads', type=int,
                        required=False, default=1,
                        dest='threads',
                        help='Number of threads used for download. You should '
                             'provide an API key to perform more requests '
                             'through Entrez.')

    parser.add_argument('-r', '--retry', type=int,
                        required=False, dest='retry',
                        default=7,
                        help='Maximum number of retries when a '
                             'download fails.')

    parser.add_argument('-k', '--api-key', type=str,
                        required=False, dest='api_key',
                        help='Personal API key provided to the NCBI. If not set, '
                             'only 3 requests per second are allowed. With a '
                             'valid API key the limit increases to 10 '
                             'requests per second.')

    parser.add_argument('-fm', '--fetch-metadata',
                        required=False, dest='fetch_metadata',
                        action='store_true',
                        default=False,
                        help='If provided, the process downloads '
                             'metadata for the assemblies.')

    parser.add_argument('-f', '--filtering-criteria', type=pv.validate_criteria_file,
                        required=False, dest='filtering_criteria',
                        help='TSV file containing filtering parameters '
                             'applied before assembly download.')

    parser.add_argument('--download', action='store_true',
                        required=False, dest='download',
                        help='If the assemblies that passed the filtering '
                             'criteria should be downloaded.')

    # Arguments specific for NCBI
    parser.add_argument('-i', '--input-table', type=str,
                        required=False, dest='input_table',
                        help='Text file with a list of accession numbers for the NCBI Assembly database.')

    args = parser.parse_args()

    del args.DownloadAssemblies

    DownloadAssemblies.main(args)

def refine_schema():
    
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    
    parser.add_argument('RefineSchema', nargs='+',
                        help='')

    parser.add_argument('-s', '--schema', type=str,
                        required=True, dest='schema',
                        help='')

    parser.add_argument('-o', '--output-directory', type=str,
                        required=True, dest='output_directory',
                        help='Path to the directory to which '
                             'files will be stored.')
    
    parser.add_argument('-a', '--allelecall-directory', type=str,
                        required=True, dest='allelecall_directory',
                        help='Path to the directory that contains'
                             'allele call directory that was run'
                             'with --no-cleanup.')
    
    parser.add_argument('-ap', '--alignment_ratio_threshold_paralagous', type=float,
                        required=False, dest='alignment_ratio_threshold_paralagous',
                        default=0.6, help='Threshold value for alignment used to '
                        'identify paralagous loci (float: 0-1).')
    
    parser.add_argument('-pp', '--pident_threshold_paralagous', type=int,
                    required=False, dest='pident_threshold_paralagous',
                    default=70, help='Threshold value for pident values used to '
                    'identify paralagous loci (int 0-100).')
    
    parser.add_argument('-ag', '--alignment_ratio_threshold_gene_fusions', type=float,
                        required=False, dest='alignment_ratio_threshold_gene_fusions',
                        default=0.9, help='Threshold value for alignment used to '
                        'indentify gene fusions (float: 0-1).')
    
    parser.add_argument('-pg', '--pident_threshold_gene_fusions', type=int,
                    required=False, dest='pident_threshold_gene_fusions',
                    default=90, help='Threshold value for pident values used to '
                    'indentify gene fusions (int 0-100).')
    
    parser.add_argument('-cs', '--clustering-sim', type=int,
                    required=False, dest='clustering_sim',
                    default=0.9, help='Similiriaty value for'
                    'kmers representatives (float: 0-1).')
    
    parser.add_argument('-cc', '--clustering-cov', type=int,
                    required=False, dest='clustering_cov',
                    default=0.9, help='Coverage value for'
                    'kmers representatives (float: 0-1).')
    
    parser.add_argument('-gp', '--genome_presence', type=int,
                    required=False, dest='genome_presence',
                    help='The minimum number of genomes specific cluster'
                    'cluster of CDS must be present in order to be considered')
    
    parser.add_argument('-c', '--cpu', type=int,
                    required=False, dest='cpu',
                    default=1, 
                    help='Number of cpus to run blast instances.')

    args = parser.parse_args()

    del args.RefineSchema

    SpuriousLociIdentification.main(**vars(args))

def main():

    module_info = {"DownloadAssemblies": ['Downloads assemblies from the NCBI '
                                       'and the ENA661K database.', download_assemblies],
                   "RefineSchema": ['Identifies spurious loci from a schema', refine_schema]}

    if len(sys.argv) == 1 or sys.argv[1] not in module_info:
        print('USAGE: SchemaRefinery [module] -h \n')
        print('Select one of the following modules:\n')
        for f in module_info:
            print('{0}: {1}'.format(f, module_info[f][0]))
        sys.exit(0)

    module = sys.argv[1]
    module_info[module][1]()


if __name__ == "__main__":

    main()
