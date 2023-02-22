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
except ModuleNotFoundError:
    from Schema_refinery.DownloadAssemblies import download_assemblies

def download_module():
    """
    Function that contains the required arguments for the download_module of the
    schema refinery package, this function exports the arguments to
    download_asseblies.py
    """
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument('Download_module', nargs='+',
                        help='')

    #Common arguments between databases
    parser.add_argument('-db', '--database', type=str,
                        required=True, dest='database',
                        nargs='+',
                        choices = ['NCBI','ENA661K'],
                        help='Databases from which assemblies will be downloaded.')

    parser.add_argument('-o', '--output-directory', type=str,
                        required=True, dest='output_directory',
                        help='Path to the output directory.')

    parser.add_argument('-e', '--email', type=str,
                        required=True, dest='email',
                        help='Entrez email parameter.')

    parser.add_argument('-t', '--taxon', type=str,
                        required=False, dest='taxon',
                        help='Scientific name of the taxon.')

    parser.add_argument('-th', '--threads', type=int,
                        required=False, default=3,
                        dest='threads',
                        help='Number of threads for download.')

    parser.add_argument('-r', '--retry', type=int,
                        required=False, dest='retry',
                        default=7,
                        help='Maximum number of retries when a '
                             'download fails.')

    parser.add_argument('-k', '--api-key', type=str,
                        required=False, dest='api_key',
                        help='Personal API key from NCBI. '
                        'If not set, only 3 queries per second are allowed. '
                        '10 queries per second otherwise with a valid API key.')

    parser.add_argument('-fm', '--fetch-metadata',
                        required=False, dest='f_metadata',
                        action='store_true',
                        default = False,
                        help='Do not fetch metadata if toggled')

    parser.add_argument('-f', '--filter-criteria-path',type=str,
                        required=False, dest='filter_criteria_path',
                        help='TSV file containing filter parameters for choosen'
                             'assemblies before downloading')

    parser.add_argument('--download', action='store_true',
                        required=False, dest='download',
                        help='If the assemblies from the selected samples'
                             'should be downloaded.')

    #Arguments specific for NCBI

    parser.add_argument('-i', '--input-table', type=str,
                        required=False, dest='input_table',
                        help='TSV file downloaded from the '
                             'NCBI Genome Assembly and Annotation '
                             'report.')

    args = parser.parse_args()

    del args.Download_module

    download_assemblies.main(args)

def main():
    module_info = {"Download_module":['Downloads assemblies from either NCBI '
                                      'or ENA661K database',download_module]}

    if len(sys.argv) == 1 or sys.argv[1] not in module_info:
        print('possible arguments here')
        sys.exit()

    module = sys.argv[1]
    module_info[module][1]()



if __name__ == "__main__":

    main()
