"""Main module."""

import argparse
import sys

try:
    from __init__ import __version__
    from SchemaAnnotation import schema_annotation
    from utils import parameters_validation as pv
except:
    from schema_refinery.__init__ import __version__
    from SchemaAnnotation.SchemaAnnotation import schema_annotation
    from schema_refinery.utils import parameters_validation as pv

version = __version__

def schema_annotation_module():
    parser = argparse.ArgumentParser(prog='SchemaAnnotation',
                                     description='')

    parser.add_argument('SchemaAnnotation', nargs='+',
                        help='')

    parser.add_argument('-ao', '--annotation-options', type=str,
                        required=True, dest='annotation_options',
                        nargs='+',
                        choices = ['proteomes','genbank', 'uniprot'],
                        help='Annotation options to run '
                             'proteomes, genbank and/or uniprot')

    parser.add_argument('-ps', '--prot_species', type=str, required=False,
                        dest='species',
                        help='Uniprot Finder output for species')

    parser.add_argument('-pg', '--prot_genus', type=str, required=False,
                        dest='genus',
                        default = '',
                        help='Uniprot Finder output for genus')

    parser.add_argument('-i', type=str, required=False,
                            dest='input_files',
                            help='Path to the directory that contains '
                                'Genbank files with annotations to '
                                'extract.')

    parser.add_argument('-s', '--schema-directory', type=str,
                            required=False, dest='schema_directory',
                            help='Path to the schema\'s directory.')

    parser.add_argument('-t', '--input_table', type=str,
                        required=False, dest='input_table',
                        help='TSV file downloaded from UniProt '
                            'that contains list of proteomes.')

    parser.add_argument('-d', type=str, required=False,
                            dest='proteomes_directory',
                            help='Path to directory with UniProt '
                                'proteomes in Fasta format.')

    parser.add_argument('-th', '--threads', type=int,
                        required=False, default=2,
                        dest='threads',
                        help='Number of threads for concurrent download.')

    parser.add_argument('-r', '--retry', type=int,
                        required=False, dest='retry',
                        default=7,
                        help='Maximum number of retries when a '
                            'download fails.')
    
    parser.add_argument('-cpu', '--cpu-cores', type=int, required=False,
                            dest='cpu_cores',
                            default=1,
                            help='Number of CPU cores to pass to BLAST.')

    parser.add_argument('-o', '--output-directory', type=str,
                        required=False, dest='output_directory',
                        help='Path to the output directory.')

    args = parser.parse_args()
    del args.SchemaAnnotation

    schema_annotation(args)

def main():

    modules_info = {'SchemaAnnotation': ['This module handles the whole process of creating the'
                                        'Genbank, TrEMBL and Swiss-Prot annotations and merging them'
                                        'alongside the UniProtFinder annotations into a single TSV file.',
                                       schema_annotation_module],
    }

    matches = ["--v", "-v", "-version", "--version"]
    if len(sys.argv) > 1 and any(m in sys.argv[1] for m in matches):
        # print version and exit
        print(f'Schema Refinery version: {version}')
        sys.exit(0)

    print(f'\nSchema Refinery version: {version}')

    # display help message if selected process is not valid
    if len(sys.argv) == 1 or sys.argv[1] not in modules_info:
        print('USAGE: schema_refinery.py [module] -h \n')
        print('Select one of the following modules :\n')
        for f in modules_info:
            print(f'{f}: {modules_info[f][0]}')
        sys.exit(0)

    # Check python version
    python_version = pv.validate_python_version()

    process = sys.argv[1]
    modules_info[process][1]()


if __name__ == '__main__':

    main()
