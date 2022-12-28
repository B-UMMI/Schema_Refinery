"""Main module."""

import os
import sys
import pickle
import shutil
import hashlib
import argparse

try:
    from Download_module import Download_module

except ModuleNotFoundError:
    from Schema_refinery.Download_module import Download_module
    
def download_module():
    
    def msg(name=None):
        #Download assembles from NCBI
        usage_msg = "write_stuff"
        return usage_msg

    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument('Download_module', nargs='+',
                        help='')

    #Common arguments between databases
    parser.add_argument('-db', '--database', type=str,
                        required=True, dest='database',
                        choices = ['NCBI','ENA661K'],
                        help='Database where to fetch assemblies '
                             'ENA661K or NCBI')

    parser.add_argument('-s', '--species', type=str,
                        required=True, dest='species',
                        help='Desired species to fetch')

    parser.add_argument('-o', '--output-directory', type=str,
                        required=True, dest='output_directory',
                        help='Output Path')
    
    parser.add_argument('-th', '--threads', type=int,
                        required=False, default=2,
                        dest='threads',
                        help='Number of threads for download.')

    parser.add_argument('-r', '--retry', type=int,
                        required=False, dest='retry',
                        default=7,
                        help='Maximum number of retries when a '
                             'download fails.')
    
    parser.add_argument('-e', '--email', type=str,
                        required=True, dest='email',
                        help='email for entrez NCBI')
    
    parser.add_argument('-k', '--api_key', type=str,
                        required=False, dest='api_key',
                        help='API key to increase the mumber of requests')
    
    parser.add_argument('-fm', '--f_metadata',
                        required=False, dest='f_metadata',
                        action='store_false',
                        default = True,
                        help='Do not fetch metadata if toggled')
    
    parser.add_argument('-f', '--filter_criteria_path',type=str,
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

    #Specific for ENA661k database
    
    parser.add_argument('-p', '--paths_ena', type=str,
                        required=False, dest='paths_table',
                        help='Maximum number of retries when a '
                             'download fails.')
    
    parser.add_argument('-m', '--metadata_ena', type=str,
                        required=False, dest='metadata_table',
                        help='Maximum number of retries when a '
                             'download fails.')

    parser.add_argument('-stride', '--stride', type=str,
                        required=False,
                        dest='stride',
                        help='Interval specifying which sample ids to '
                             'download. Example: "1:2000" - This will '
                             'download the first 2000 samples. Note: If '
                             'you want to download from the first id, '
                             'you have to put "1", not "0" in the lower '
                             'value.')
    
    args = parser.parse_args()
    
    del args.Download_module

    Download_module.main(args)

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
