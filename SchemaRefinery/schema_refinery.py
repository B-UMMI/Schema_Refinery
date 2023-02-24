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
    from DownloadAssemblies import download_assemblies
    from DownloadAssemblies import parameter_validation as pv
except ModuleNotFoundError:
    from Schema_refinery.DownloadAssemblies import download_assemblies
    from Schema_refinery.DownloadAssemblies import parameter_validation as pv


def download_module():

    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument('Download_module', nargs='+',
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
                        required=False, default=3,
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
                        help='Personal API key provided to NCBI. If not set, '
                             'only 3 requests per second are allowed. With a '
                             'valid API key the limit increases to 10 '
                             'requests per second.')

    parser.add_argument('-fm', '--fetch-metadata',
                        required=False, dest='fetch_metadata',
                        action='store_true',
                        default=False,
                        help='If provided, the process does not fetch '
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
                        help='TSV file downloaded from the NCBI Genome '
                             'Assembly and Annotation report.')

    args = parser.parse_args()

    del args.Download_module

    download_assemblies.main(args)


def main():

    module_info = {"Download_module": ['Downloads assemblies from either NCBI '
                                       'or ENA661K database', download_module]}

    if len(sys.argv) == 1 or sys.argv[1] not in module_info:
        print('possible arguments here')
        sys.exit()

    module = sys.argv[1]
    module_info[module][1]()


if __name__ == "__main__":

    main()
