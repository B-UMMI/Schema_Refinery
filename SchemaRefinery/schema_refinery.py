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
    from SchemaAnnotation import SchemaAnnotation
    from utils import parameter_validation as pv
    from RefineSchema import IdentifySpuriousGenes
    from IdentifyParalagousLoci import IdentifyParalogousLoci
    from IdentifyDuplicateGenes import IdentifyDuplicateGenes
    from AdaptLoci import AdaptLoci
    from MatchSchema import MatchSchemas
    from utils import constants as ct
     
except ModuleNotFoundError:
    from SchemaRefinery.DownloadAssemblies import DownloadAssemblies
    from SchemaRefinery.SchemaAnnotation import SchemaAnnotation
    from SchemaRefinery.utils import parameter_validation as pv
    from SchemaRefinery.RefineSchema import IdentifySpuriousGenes
    from SchemaRefinery.IdentifyParalagousLoci import IdentifyParalogousLoci
    from SchemaRefinery.IdentifyDuplicateGenes import IdentifyDuplicateGenes
    from SchemaRefinery.AdaptLoci import AdaptLoci
    from SchemaRefinery.MatchSchema import MatchSchemas
    from SchemaRefinery.utils import constants as ct


def download_assemblies() -> None:
    """
    Parse command-line arguments and initiate the download of assemblies.

    This function sets up an argument parser to handle various command-line
    options for downloading assemblies from specified databases. It then
    calls the main function of the DownloadAssemblies class with the parsed
    arguments.

    Parameters
    ----------
    None

    Returns
    -------
    None
    """
    
    # Initialize the argument parser
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)

    # Common arguments between databases
    parser.add_argument('-db',
                        '--database',
                        type=str,
                        required=True,
                        dest='database',
                        nargs='+',
                        choices=ct.DATABASE_CHOICES,
                        help='Databases from which assemblies will be downloaded.')

    parser.add_argument('-o',
                        '--output-directory',
                        type=str,
                        required=True,
                        dest='output_directory',
                        help='Path to the output directory.')

    parser.add_argument('-e',
                        '--email',
                        type=str,
                        required=True,
                        dest='email',
                        help='Email provided to Entrez.')

    parser.add_argument('-t',
                        '--taxon',
                        type=str,
                        required=False,
                        dest='taxon',
                        help='Scientific name of the taxon.')

    parser.add_argument('-th',
                        '--threads',
                        type=int,
                        required=False,
                        default=1,
                        dest='threads',
                        help='Number of threads used for download. You should provide an API'
                        'key to perform more requests through Entrez.')

    parser.add_argument('-r',
                        '--retry',
                        type=int,
                        required=False,
                        dest='retry',
                        default=7,
                        help='Maximum number of retries when a download fails.')

    parser.add_argument('-k',
                        '--api-key',
                        type=str,
                        required=False,
                        dest='api_key',
                        help='Personal API key provided to the NCBI. If not set, only 3 requests per second are allowed. With a valid API key the limit increases to 10 requests per second.')

    parser.add_argument('-fm',
                        '--fetch-metadata',
                        required=False,
                        dest='fetch_metadata',
                        action='store_true',
                        default=False,
                        help='If provided, the process downloads metadata for the assemblies.')

    parser.add_argument('-f',
                        '--filtering-criteria',
                        type=pv.validate_criteria_file,
                        required=False,
                        dest='filtering_criteria',
                        help='TSV file containing filtering parameters applied before assembly download.')

    parser.add_argument('--download',
                        action='store_true',
                        required=False,
                        dest='download',
                        help='If the assemblies that passed the filtering criteria should be downloaded.')

    # Arguments specific for NCBI
    parser.add_argument('-i',
                        '--input-table',
                        type=str,
                        required=False,
                        dest='input_table',
                        help='Text file with a list of accession numbers for the NCBI Assembly database.')

    # Parse the command-line arguments
    args = parser.parse_args()

    # Call the main function of the DownloadAssemblies class with the parsed arguments
    DownloadAssemblies.main(args)


def schema_annotation() -> None:
    """
    Parse command-line arguments and initiate the schema annotation process.

    This function sets up an argument parser to handle various command-line
    options for annotating schemas. It then calls the main function of the
    SchemaAnnotation class with the parsed arguments.

    Parameters
    ----------
    None

    Returns
    -------
    None
    """
    
    # Initialize the argument parser
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)

    # Add arguments to the parser
    parser.add_argument('-s',
                        '--schema-directory',
                        type=str,
                        required=True,
                        dest='schema_directory',
                        help='Path to the schema\'s directory.')

    parser.add_argument('-o',
                        '--output-directory',
                        type=str,
                        required=True,
                        dest='output_directory',
                        help='Path to the output directory where to save the files.')

    parser.add_argument('-ao',
                        '--annotation-options',
                        type=str,
                        required=True,
                        dest='annotation_options',
                        nargs='+',
                        choices=ct.SCHEMA_ANNOTATION_RUNS_CHOICES,
                        help='Annotation options to run. "uniprot-proteomes" to download UniProt reference proteomes for the taxa and align with BLASTp. "genbank-files" to align against the CDSs in a set of Genbank files. "uniprot-sparql" to search for exact matches through UniProt\'s SPARQL endpoint. "match-schemas" to align against provided target schema and report best matches.')

    parser.add_argument('-pt',
                        '--proteome-table',
                        type=str,
                        required=False,
                        dest='proteome_table',
                        help='TSV file downloaded from UniProt that contains the list of proteomes.')

    parser.add_argument('-gf',
                        '--genbank-files',
                        type=str,
                        required=False,
                        dest='genbank_files',
                        help='Path to the directory that contains Genbank files with annotations to extract.')

    parser.add_argument('-ca',
                        '--chewie-annotations',
                        type=str,
                        required=False,
                        nargs='+',
                        dest='chewie_annotations',
                        help='File with the results from chewBBACA UniprotFinder module.')

    parser.add_argument('-ss',
                        '--subject-schema',
                        type=str,
                        required=False,
                        dest='subject_schema',
                        help='Path to the subject schema directory. This argument is needed by the Match Schemas sub-module.')

    parser.add_argument('-sa',
                        '--subject-annotations',
                        type=str,
                        required=False,
                        dest='subject_annotations',
                        help='Annotations of the subject schema.')

    parser.add_argument('-sc',
                        '--subject-columns',
                        type=str,
                        required=False,
                        nargs='+',
                        dest='subject_columns',
                        help='Columns from the subject schema annotations to merge into the new schema annotations.')

    parser.add_argument('--bsr',
                        type=float,
                        required=False,
                        default=0.6,
                        dest='bsr',
                        help='Minimum BSR value to consider aligned alleles as alleles for the same locus. This argument is optional for the Match Schemas sub-module.')
    
    parser.add_argument('-t',
                        '--threads',
                        type=int,
                        required=False,
                        default=1,
                        dest='threads',
                        help='Number of threads for concurrent download.')

    parser.add_argument('-c',
                        '--cpu',
                        type=int,
                        required=False,
                        default=1,
                        dest='cpu',
                        help='Number of CPU cores for multiprocessing.')

    parser.add_argument('-r',
                        '--retry',
                        type=int,
                        required=False,
                        dest='retry',
                        default=7,
                        help='Maximum number of retries when a download fails.')
    
    parser.add_argument('-tt',
                        '--translation-table',
                        type=int,
                        required=False,
                        dest='translation_table',
                        default=11,
                        help='Translation table to use for the CDS translation.')

    parser.add_argument('-cs',
                        '--clustering-sim',
                        type=float,
                        required=False,
                        dest='clustering_sim',
                        default=0.9,
                        help='Similarity value for kmers representatives (float: 0-1).')

    parser.add_argument('-cc',
                        '--clustering-cov',
                        type=float,
                        required=False,
                        dest='clustering_cov',
                        default=0.9,
                        help='Coverage value for kmers representatives (float: 0-1).')
    
    parser.add_argument('-sr',
                        '--size_ratio',
                        type=float,
                        required=False,
                        dest='size_ratio',
                        default=0.8,
                        help='Size ratio to consider alleles as the same locus.')
    
    parser.add_argument('-rm',
                        '--run-mode',
                        type=str,
                        required=False,
                        dest='run_mode',
                        default='reps',
                        choices=ct.SCHEMA_ANNOTATION_RUN_MODE_CHOICES,
                        help='Mode to run the module: reps or alleles.')

    parser.add_argument('-pm',
                        '--processing-mode',
                        type=str,
                        required=False,
                        dest='processing_mode',
                        default='reps_vs_alleles',
                        choices=ct.PROCESSING_MODE_CHOICES,
                        help='Mode to run the module for Schema match: reps_vs_reps,'
                        'reps_vs_alleles, alleles_vs_alleles, alleles_vs_reps.')
    
    parser.add_argument('-egtc',
                        '--extra_genbank_table_columns',
                        type=str,
                        required=False,
                        dest='extra_genbank_table_columns',
                        nargs='+',
                        default=[],
                        choices=ct.GENBANK_CDS_QUALIFIERS_CHOICES,
                        help='List of columns to add to annotation file (locus_tag, note, codon_start, function, protein_id, db_xref).')

    parser.add_argument('-gia',
                    '--genbank-ids-to-add',
                    type=str,
                    required=False,
                    dest='genbank_ids_to_add',
                    nargs='+',
                    default=[],
                    choices=ct.GENBANK_CDS_QUALIFIERS_CHOICES,
                    help='List of GenBank IDs to add to final results.')
    
    parser.add_argument('-pia',
                    '--proteome-ids-to-add',
                    type=str,
                    required=False,
                    dest='proteome_ids_to_add',
                    nargs='+',
                    default=[],
                    choices=ct.GENBANK_CDS_QUALIFIERS_CHOICES,
                    help='List of Proteome IDs to add to final results.')

    parser.add_argument('--nocleanup',
                        action='store_true',
                        required=False,
                        dest='no_cleanup',
                        help='Flag to indicate whether to skip cleanup after running the module.')

    # Parse the command-line arguments
    args = parser.parse_args()

    # Call the main function of the SchemaAnnotation class with the parsed arguments
    SchemaAnnotation.main(args)


def identify_spurious_genes() -> None:
    """
    Parse command-line arguments and initiate the process to identify spurious genes.

    This function sets up an argument parser to handle various command-line
    options for identifying spurious genes in a schema. It then calls the main
    function of the IdentifySpuriousGenes class with the parsed arguments.

    Parameters
    ----------
    None

    Returns
    -------
    None
    """
    
    # Initialize the argument parser
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)

    # Add arguments to the parser
    parser.add_argument('-s',
                        '--schema-directory',
                        type=str,
                        required=True,
                        dest='schema_directory',
                        help='Path to the created schema directory.')

    parser.add_argument('-o',
                        '--output-directory',
                        type=str,
                        required=True,
                        dest='output_directory',
                        help='Path to the directory to which files will be stored.')

    parser.add_argument('-a',
                        '--allelecall-directory',
                        type=str,
                        required=True,
                        dest='allelecall_directory',
                        help='Path to the directory that contains allele call directory that was run with --no-cleanup.')
    
    parser.add_argument('-pnl',
                        '--possible-new-loci',
                        type=str,
                        required=False,
                        dest='possible_new_loci',
                        help='Path to the directory that contains possible new loci.')

    parser.add_argument('-at',
                        '--alignment_ratio_threshold',
                        type=float,
                        required=False,
                        dest='alignment_ratio_threshold',
                        default=0.9,
                        help='Threshold value for alignment used to identify spurious CDS (float: 0-1).')

    parser.add_argument('-pt',
                        '--pident_threshold',
                        type=int,
                        required=False,
                        dest='pident_threshold',
                        default=90,
                        help='Threshold value for pident values used to identify spurious CDS (int 0-100).')

    parser.add_argument('-cs',
                        '--clustering-sim',
                        type=float,
                        required=False,
                        dest='clustering_sim',
                        default=0.9,
                        help='Similarity value for kmers representatives (float: 0-1).')

    parser.add_argument('-cc',
                        '--clustering-cov',
                        type=float,
                        required=False,
                        dest='clustering_cov',
                        default=0.9,
                        help='Coverage value for kmers representatives (float: 0-1).')

    parser.add_argument('-gp',
                        '--genome_presence',
                        type=int,
                        required=False,
                        dest='genome_presence',
                        help='The minimum number of genomes specific cluster of CDS must be present in order to be considered.')

    parser.add_argument('-as',
                        '--absolute_size',
                        type=int,
                        required=False,
                        dest='absolute_size',
                        default=201,
                        help='Size of the CDS to consider processing.')

    parser.add_argument('-tt',
                        '--translation_table',
                        type=int,
                        required=False,
                        dest='translation_table',
                        default=11,
                        help='Translation table to use for the CDS translation.')
    
    parser.add_argument('-b',
                        '--bsr',
                        type=float,
                        required=False,
                        dest='bsr',
                        default=0.6,
                        help='BSR value to consider alleles as the same locus.')
    
    parser.add_argument('-sr',
                        '--size_ratio',
                        type=float,
                        required=False,
                        dest='size_ratio',
                        default=0.8,
                        help='Size ratio to consider alleles as the same locus.')
    
    parser.add_argument('-m',
                        '--run-mode',
                        type=str,
                        required=False,
                        dest='run_mode',
                        default='schema',
                        choices=ct.IDENTIFY_SPURIOUS_LOCI_RUN_MODE_CHOICES,
                        help='Run mode for identifying spurious loci.')

    parser.add_argument('-pm',
                        '--processing-mode',
                        type=str,
                        required=False,
                        dest='processing_mode',
                        default='reps_vs_alleles',
                        choices=ct.PROCESSING_MODE_CHOICES,
                        help='Mode to run the module: reps_vs_reps, reps_vs_alleles, alleles_vs_alleles, alleles_vs_reps.')

    parser.add_argument('-c',
                        '--cpu',
                        type=int,
                        required=False,
                        dest='cpu',
                        default=1, 
                        help='Number of CPUs to run BLAST instances.')

    parser.add_argument('--nocleanup',
                        action='store_true',
                        required=False,
                        dest='no_cleanup',
                        help='Flag to indicate whether to skip cleanup after running the module.')

    # Parse the command-line arguments
    args = parser.parse_args()

    # Validate the possible_new_loci argument
    if args.possible_new_loci and args.run_mode != 'loci_vs_cds':
        sys.exit("Argument -p --possible_new_loci can only be used with -m --run-mode of loci_vs_cds.")
    
    # Call the main function of the IdentifySpuriousGenes class with the parsed arguments
    IdentifySpuriousGenes.main(**vars(args))


def adapt_loci() -> None:
    """
    Parse command-line arguments and initiate the process to adapt loci.

    This function sets up an argument parser to handle various command-line
    options for adapting loci. It then calls the main function of the AdaptLoci
    class with the parsed arguments.

    Parameters
    ----------
    None

    Returns
    -------
    None
    """
    
    # Initialize the argument parser
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)

    # Add arguments to the parser
    parser.add_argument('-i',
                        '--input_file',
                        type=str,
                        required=True,
                        dest='input_file',
                        help='TSV file with the loci path to be adapted.')
    
    parser.add_argument('-o',
                        '--output-directory',
                        type=str,
                        required=True,
                        dest='output_directory',
                        help='Path to the directory to which files will be stored.')
    
    parser.add_argument('-c',
                        '--cpu',
                        type=int,
                        required=False,
                        dest='cpu_cores',
                        default=1, 
                        help='Number of CPUs to run BLAST instances.')
    
    parser.add_argument('-b',
                        '--bsr',
                        type=float,
                        required=False,
                        dest='blast_score_ratio',
                        default=0.6,
                        help='BSR value to consider alleles as the same locus.')
    
    parser.add_argument('-tt',
                        '--translation_table',
                        type=int,
                        required=False,
                        dest='translation_table',
                        default=11,
                        help='Translation table to use for the CDS translation.')
    
    # Parse the command-line arguments
    args = parser.parse_args()

    # Call the main function of the AdaptLoci class with the parsed arguments
    AdaptLoci.main(**vars(args))


def identify_paralogous_loci() -> None:
    """
    Parse command-line arguments and initiate the process to identify paralogous loci.

    This function sets up an argument parser to handle various command-line
    options for identifying paralogous loci in a schema. It then calls the main
    function of the IdentifyParalogousLoci class with the parsed arguments.

    Parameters
    ----------
    None

    Returns
    -------
    None
    """
    
    # Initialize the argument parser
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    
    # Add arguments to the parser
    parser.add_argument('-s',
                        '--schema-directory',
                        type=str,
                        required=True,
                        dest='schema_directory',
                        help='Folder that contains the schema to identify paralogous loci.')
    
    parser.add_argument('-o',
                        '--output-directory',
                        type=str,
                        required=True,
                        dest='output_directory',
                        help='Path to the directory to which files will be stored.')
    
    parser.add_argument('-c',
                        '--cpu',
                        type=int,
                        required=False,
                        dest='cpu',
                        default=1, 
                        help='Number of CPUs to run BLAST instances.')
    
    parser.add_argument('-b',
                        '--bsr',
                        type=float,
                        required=False,
                        dest='bsr',
                        default=0.6,
                        help='BSR value to consider alleles as the same locus.')
    
    parser.add_argument('-tt',
                        '--translation_table',
                        type=int,
                        required=False,
                        dest='translation_table',
                        default=11,
                        help='Translation table to use for the CDS translation.')
    
    parser.add_argument('-st',
                        '--size_threshold',
                        type=float,
                        required=False,
                        dest='size_threshold',
                        default=0.2,
                        help="Size threshold to consider two paralogous loci as similar.")
    
    parser.add_argument('-pm',
                        '--processing-mode',
                        type=str,
                        required=False,
                        dest='processing_mode',
                        choices=ct.PROCESSING_MODE_CHOICES,
                        default='alleles_vs_alleles',
                        help='Mode to run the module: reps_vs_reps, reps_vs_alleles, alleles_vs_alleles, alleles_vs_reps.')

    parser.add_argument('--nocleanup',
                        action='store_true',
                        required=False,
                        dest='no_cleanup',
                        help='Flag to indicate whether to skip cleanup after running the module.')

    # Parse the command-line arguments
    args = parser.parse_args()

    # Call the main function of the IdentifyParalogousLoci class with the parsed arguments
    IdentifyParalogousLoci.identify_paralogous_loci(**vars(args))


def identify_duplicate_gene() -> None:
    """
    Parse command-line arguments and initiate the process to identify problematic loci.

    This function sets up an argument parser to handle various command-line
    options for identifying problematic loci in a schema. It then calls the main
    function of the IdentifyProblematicLoci class with the parsed arguments.

    Parameters
    ----------
    None

    Returns
    -------
    None
    """
    
    # Initialize the argument parser
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    
    # Add arguments to the parser
    parser.add_argument('-d',
                        '--distinct-hashtable',
                        type=str,
                        required=True,
                        dest='distinct_hashtable',
                        help='Path to the distinct.hashtable file containing CDS sequences.')
    
    parser.add_argument('-s',
                        '--schema-directory',
                        type=str,
                        required=True,
                        dest='schema_directory',
                        help='Path to the directory containing schema loci.')
    
    parser.add_argument('-o',
                        '--output-directory',
                        type=str,
                        required=True,
                        dest='output_directory',
                        help='Path to the directory where the output file will be saved.')
    
    parser.add_argument('-pt',
                        '--problematic-threshold',
                        type=float,
                        required=False,
                        dest='problematic_threshold',
                        default=0.1,
                        help='Threshold for determining if a locus is problematic.')

    parser.add_argument('--nocleanup',
                        action='store_true',
                        required=False,
                        dest='no_cleanup',
                        help='Flag to indicate whether to skip cleanup after running the module.')

    # Parse the command-line arguments
    args = parser.parse_args()

    # Call the main function of the IdentifyProblematicLoci class with the parsed arguments
    IdentifyDuplicateGenes.identify_duplicate_gene(**vars(args))


def match_schemas() -> None:
    """
    Parse command-line arguments and initiate the process to match schemas.

    This function sets up an argument parser to handle various command-line
    options for matching schemas. It then calls the main function of the
    MatchSchemas class with the parsed arguments.

    Parameters
    ----------
    None

    Returns
    -------
    None
    """
    
    # Initialize the argument parser
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
     
    # Add arguments to the parser
    parser.add_argument('-qs',
                        '--query-schema-directory',
                        type=str,
                        required=True,
                        dest='query_schema_directory',
                        help='Folder that contains the schema to identify paralogous loci.')

    parser.add_argument('-ss',
                        '--subject-chema-directory',
                        type=str,
                        required=True,
                        dest='subject_schema_directory',
                        help='Folder that contains the schema to identify paralogous loci.')
    
    parser.add_argument('-o',
                        '--output-directory',
                        type=str,
                        required=True,
                        dest='output_directory',
                        help='Path to the directory to which files will be stored.')
    
    parser.add_argument('-c',
                        '--cpu',
                        type=int,
                        required=False,
                        dest='cpu',
                        default=1, 
                        help='Number of CPUs to run BLAST instances.')
    
    parser.add_argument('-b',
                        '--bsr',
                        type=float,
                        required=False,
                        dest='bsr',
                        default=0.6,
                        help='BSR value to consider alleles as the same locus.')
    
    parser.add_argument('-tt',
                        '--translation_table',
                        type=int,
                        required=False,
                        dest='translation_table',
                        default=11,
                        help='Translation table to use for the CDS translation.')

    parser.add_argument('-pm',
                        '--processing-mode',
                        type=str,
                        required=False,
                        dest='processing_mode',
                        choices=ct.PROCESSING_MODE_CHOICES,
                        default='alleles_vs_alleles',
                        help='Mode to run the module: reps_vs_reps, reps_vs_alleles, alleles_vs_alleles, alleles_vs_reps.')
    
    parser.add_argument('--nocleanup',
                        action='store_true',
                        required=False,
                        dest='no_cleanup',
                        help='Flag to indicate whether to skip cleanup after running the module.')

    # Parse the command-line arguments
    args = parser.parse_args()

    # Call the main function of the MatchSchemas class with the parsed arguments
    MatchSchemas.match_schemas(**vars(args))
    

def main():

    module_info = {'DownloadAssemblies': ["Downloads assemblies from the NCBI "
                                    "and the ENA661K database.", download_assemblies],
                        'SchemaAnnotation': ['Annotate a schema based on TrEMBL and Swiss-Prot '
                                            'records, and based on alignment against Genbank '
                                            'files and other schemas.',
                                            schema_annotation],
                        'IdentifySpuriousGenes': ["Identifies spurious genes in a schema by running against itself or"
                                        " against unclassified CDS to infer new loci and identify problematic genes.",
                                        identify_spurious_genes],
                        'AdaptLoci': ["Adapts loci from a fasta files to a new schema.", adapt_loci],
                        'IdentifyParalagousLoci': ["Identifies paralagous loci based on schema input", identify_paralogous_loci],
                        'IdentifyDuplicateGenes': ["Identify problematic loci based on the presence of NIPHs and NIPHEMs.", identify_duplicate_gene],
                        'MatchSchema': ["Match schemas to identify the best matches between two schemas.", match_schemas]}

    if len(sys.argv) == 1 or sys.argv[1] not in module_info:
        print('USAGE: SchemaRefinery [module] -h \n')
        print('Select one of the following modules:\n')
        for f in module_info:
            print('{0}: {1}'.format(f, module_info[f][0]))
        sys.exit(0)

    module = sys.argv[1]
    sys.argv.remove(module)
    module_info[module][1]()

if __name__ == "__main__":

    main()
