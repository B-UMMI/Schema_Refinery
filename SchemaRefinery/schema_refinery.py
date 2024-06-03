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

from RefineSchema import SpuriousLoci

try:
    from DownloadAssemblies import DownloadAssemblies
    from SchemaAnnotation import SchemaAnnotation
    from utils import parameter_validation as pv
    from RefineSchema import (UnclassifiedCDS,
                              SpuriousLoci)
except ModuleNotFoundError:
    from SchemaRefinery.DownloadAssemblies import DownloadAssemblies
    from SchemaRefinery.SchemaAnnotation import SchemaAnnotation
    from SchemaRefinery.utils import parameter_validation as pv
    from SchemaRefinery.RefineSchema import (UnclassifiedCDS,
                                            SpuriousLoci)


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

def schema_annotation():

    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument('SchemaAnnotation', nargs='+',
                        help='')

    parser.add_argument('-s', '--schema-directory', type=str,
                        required=True, dest='schema_directory',
                        help='Path to the schema\'s directory.')

    parser.add_argument('-o', '--output-directory', type=str,
                        required=True, dest='output_directory',
                        help='Path to the output directory where to '
                             'save the files.')

    parser.add_argument('-ao', '--annotation-options', type=str,
                        required=True, dest='annotation_options',
                        nargs='+',
                        choices=['uniprot-proteomes', 'genbank',
                                 'uniprot-sparql', 'match-schemas'],
                        help='Annotation options to run. "uniprot-proteomes" '
                             'to download UniProt reference proteomes for '
                             'the taxa and align with BLASTp. "genbank-files"'
                             ' to aligned against the CDSs in a set of '
                             'Genbank files. "uniprot-sparql" to search for '
                             'exact matches through UniProt\'s SPARQL '
                             'endpoint. "match-schemas" to align against '
                             'provided target schema and report best matches.')

    parser.add_argument('-pt', '--proteome-table', type=str,
                        required=False, dest='proteome_table',
                        help='TSV file downloaded from UniProt '
                             'that contains the list of proteomes.')

    parser.add_argument('-gf', '--genbank-files', type=str,
                        required=False, dest='genbank_files',
                        help='Path to the directory that contains '
                             'Genbank files with annotations to '
                             'extract.')

    parser.add_argument('-ca', '--chewie-annotations', type=str,
                        required=False, nargs='+', dest='chewie_annotations',
                        help='File with the results from chewBBACA '
                             'UniprotFinder module.')

    parser.add_argument('-ss', '--subject-schema', type=str,
                        required=False, dest='subject_schema',
                        help='Path to que subject schema directory. '
                             'This argument is needed by the Match Schemas '
                             'sub-module.')

    parser.add_argument('-sa', '--subject-annotations', type=str,
                        required=False, dest='subject_annotations',
                        help='Annotations of the subject schema.')

    parser.add_argument('-sc', '--subject-columns', type=str, required=False,
                        nargs='+',
                        dest='subject_columns',
                        help='Columns from the subject schema annotations '
                             'to merge into the new schema annotations.')

    parser.add_argument('--bsr', type=float, required=False,
                        default=0.6, dest='blast_score_ratio',
                        help='Minimum BSR value to consider aligned '
                             'alleles as alleles for the same locus. '
                             'This argument is optional for the Match Schemas '
                             'sub-module.')

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
                        dest='cpu_cores', default=1,
                        help='Number of CPU cores to pass to BLAST.')

    args = parser.parse_args()

    del args.SchemaAnnotation

    SchemaAnnotation.main(args)

def unclassified_cds():
    
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    
    parser.add_argument('RefineSchema', nargs='+',
                        help='')

    parser.add_argument('-s', '--schema', type=str,
                        required=True, dest='schema',
                        help='Path to the created schema folder.')

    parser.add_argument('-o', '--output-directory', type=str,
                        required=True, dest='output_directory',
                        help='Path to the directory to which '
                             'files will be stored.')
    
    parser.add_argument('-a', '--allelecall-directory', type=str,
                        required=True, dest='allelecall_directory',
                        help='Path to the directory that contains'
                             'allele call directory that was run'
                             'with --no-cleanup.')
    
    parser.add_argument('-at', '--alignment_ratio_threshold', type=float,
                        required=False, dest='alignment_ratio_threshold_gene_fusions',
                        default=0.9, help='Threshold value for alignment used to '
                        'indentify gene fusions (float: 0-1).')
    
    parser.add_argument('-pt', '--pident_threshold', type=int,
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
                    'cluster of CDS must be present in order to be considered.')
    
    parser.add_argument('-st', '--size_threshold', type=int,
                    required=False, dest='size_threshold',
                    help='Size of the CDS to consider processing.')
    
    parser.add_argument('-c', '--cpu', type=int,
                    required=False, dest='cpu',
                    default=1, 
                    help='Number of cpus to run blast instances.')

    args = parser.parse_args()

    del args.RefineSchema

    UnclassifiedCDS.main(**vars(args))

def spurious_loci():
    
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    
    parser.add_argument('RefineSchema', nargs='+',
                        help='')

    parser.add_argument('-s', '--schema', type=str,
                        required=True, dest='schema',
                        help='Path to the created schema folder.')

    parser.add_argument('-o', '--output-directory', type=str,
                        required=True, dest='output_directory',
                        help='Path to the directory to which '
                             'files will be stored.')
    
    parser.add_argument('-a', '--allelecall-directory', type=str,
                        required=True, dest='allelecall_directory',
                        help='Path to the directory that contains'
                             'allele call directory that was run'
                             'with --no-cleanup.')
    
    parser.add_argument('-at', '--alignment_ratio_threshold', type=float,
                        required=False, dest='alignment_ratio_threshold_gene_fusions',
                        default=0.9, help='Threshold value for alignment used to '
                        'indentify gene fusions (float: 0-1).')
    
    parser.add_argument('-pt', '--pident_threshold', type=int,
                    required=False, dest='pident_threshold_gene_fusions',
                    default=90, help='Threshold value for pident values used to '
                    'indentify gene fusions (int 0-100).')
    
    parser.add_argument('-st', '--size_threshold', type=int,
                    required=False, dest='size_threshold',
                    help='Size of the CDS to consider processing.')
    
    parser.add_argument('-c', '--cpu', type=int,
                    required=False, dest='cpu',
                    default=1, 
                    help='Number of cpus to run blast instances.')

    args = parser.parse_args()

    del args.RefineSchema

    SpuriousLoci.main(**vars(args))

def main():

    module_info = {'DownloadAssemblies': ["Downloads assemblies from the NCBI "
                                       "and the ENA661K database.", download_assemblies],
                    'SchemaAnnotation': ['Annotate a schema based on TrEMBL and Swiss-Prot '
                                         'records, and based on alignment against Genbank '
                                         'files and other schemas.',
                                        schema_annotation],
                    'SpuriousLoci': ["Identifies spurious loci in a schema.", spurious_loci],
                    'UnclassifiedCDS': ["Classifies unclassified and"
                                        "missed classes CDS from a schema.", unclassified_cds]}

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
