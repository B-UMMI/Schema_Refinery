"""Main module."""

import argparse
import sys

try:
    from SchemaAnnotation import schema_annotation
except:
    from SchemaAnnotation.SchemaAnnotation import schema_annotation

def schema_annotation_module():
    parser = argparse.ArgumentParser(prog='SchemaAnnotation',
                                     description='')

    parser.add_argument('SchemaAnnotation', nargs='+',
                        help='')

    parser.add_argument('-ao', '--annotation-options', type=str,
                        required=True, dest='annotation_options',
                        nargs='+',
                        choices = ['proteomes','genbank'],
                        help='Annotation options to run '
                             'proteomes or genbank')

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

    module_info = {'SchemaAnnotation': ['',
                                       schema_annotation_module],
    }

    process = sys.argv[1]
    module_info[process][1]()


if __name__ == '__main__':

    main()
