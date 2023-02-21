"""Main module."""

import argparse
import sys

try:
    from __init__ import __version__
    from SchemaAnnotation.SchemaAnnotation import schema_annotation
    from utils.validation import validate_python_version
except:
    from Schema_refinery.__init__ import __version__
    from Schema_refinery.SchemaAnnotation.SchemaAnnotation import schema_annotation
    from Schema_refinery.utils.validation import validate_python_version

version = __version__

def schema_annotation_module():
    parser = argparse.ArgumentParser(prog='SchemaAnnotation',
                                     description='')

    parser.add_argument('SchemaAnnotation', nargs='+',
                        help='')

    parser.add_argument('-ao', '--annotation-options', type=str,
                        required=True, dest='annotation_options',
                        nargs='+',
                        choices = ['proteomes','genbank', 'uniprot', 'matchSchemas'],
                        help='Annotation options to run '
                             'proteomes, genbank, uniprot and/or matchSchemas')

    parser.add_argument('-us', '--uniprot_species', type=str, required=False,
                        dest='uniprot_species',
                        help='Uniprot Finder output file for species')

    parser.add_argument('-ug', '--uniprot_genus', type=str, required=False,
                        dest='uniprot_genus',
                        default = '',
                        help='Uniprot Finder output file for genus')

    parser.add_argument('-i', '--input-files', type=str, required=False,
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

    parser.add_argument('-d', '--proteomes-directory', type=str, required=False,
                            dest='proteomes_directory',
                            help='Path to the directory with UniProt '
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
                        required=True, dest='output_directory',
                        help='Path to the output directory where to save the files.')
    
    parser.add_argument('-qs', '--query-schema', type=str, required=False,
                        dest='query_schema',
                        help='Path to the query schema directory.'
                             'This schema will be matched against '
                             'the subject schema.'
                             'This argument is needed by the Match Schemas'
                             'sub-module.')

    parser.add_argument('-ss', '--subject-schema', type=str, required=False,
                        dest='subject_schema',
                        help='Path to que subject schema directory. '
                             'This argument is needed by the Match Schemas '
                             'sub-module.')
    
    parser.add_argument('-oc', '--old-schema-columns', type=str, required=False,
                        nargs='+',
                        dest='old_schema_columns',
                        help='Columns from the old schema annotations to merge into '
                             'the new schema annotations table being created. '
                             'This argument is needed by the Match Schemas '
                             'sub-module.')
    
    parser.add_argument('-ma', '--match_to_add', type=str, required=False,
                        dest='match_to_add',
                        default = '',
                        help='Annotation of another schema, needed with '
                             '--old-schema-columns. '
                             'This argument is needed by the Match Schemas '
                             'sub-module.')
    
    parser.add_argument('--bsr', type=float, required=False,
                        default=0.6, dest='blast_score_ratio',
                        help='Minimum BSR value to consider aligned '
                             'alleles as alleles for the same locus. '
                             'This argument is optional for the Match Schemas '
                             'sub-module.')

    args = parser.parse_args()
    del args.SchemaAnnotation

    schema_annotation(args)

def main():

    modules_info = {'SchemaAnnotation': ['This module handles the whole process of creating the '
                                        'Genbank, TrEMBL and Swiss-Prot annotations and merging them '
                                        'alongside the UniProtFinder annotations into a single TSV file. '
                                        'It also supports a matchSchemas sub-module that gets a user '
                                        'selected list of annotations from an old schema',
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
    python_version = validate_python_version()

    process = sys.argv[1]
    modules_info[process][1]()

if __name__ == '__main__':

    main()
