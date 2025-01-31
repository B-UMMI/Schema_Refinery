#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Purpose
-------

This module contains default values for Schema_refinery's
parameters.

Code documentation
------------------
"""
from typing import Tuple

# GitHub repository and contacts
REPOSITORY = 'https://github.com/B-UMMI/Schema_Refinery'
CONTACTS = 'imm-bioinfo@medicina.ulisboa.pt'

# Schema Refinery's version
VERSION = '0.3.2.4'

# Schema Refinery's Dependencies
DEPENDENCIES = [
        "numpy",
        "scipy",
        "biopython",
        "plotly",
        "requests",
        "pandas",
        "psutil",
    ]

# Dependencies version

DEPENDENCIES_VERSION = [
    "numpy >= 1.24.3, <2.0.0",
    "scipy >= 1.10.1",
    "biopython >= 1.79",
    "plotly >= 5.8.0",
    "requests >= 2.27.1",
    "pandas >= 1.5.1",
    "psutil >= 5.1.1",
    "tqdm >= 4.62.0",
    "networkx >= 2.6.0, <3.0.0"
]

# minimum Python version
MIN_PYTHON = [(3, 6, 0), '3.6.0']

# socket timeout for urllib calls
SOCKET_TIMEOUT = 30

# URL template for proteome download
PROTEOME_TEMPLATE_URL = 'https://rest.uniprot.org/uniprotkb/stream?query=proteome:{0}&format=fasta&compressed=true'

# DNA bases
DNA_BASES = ['A', 'T', 'C', 'G']

# Schema Refinery's default values

PROCESSING_MODE_CHOICES = ['reps_vs_reps', 'reps_vs_alleles', 'alleles_vs_alleles', 'alleles_vs_reps']

IDENTIFY_SPURIOUS_LOCI_RUN_MODE_CHOICES = ['unclassified_cds', 'schema']

SCHEMA_ANNOTATION_RUN_MODE_CHOICES = ['reps', 'alleles']

DATABASE_CHOICES = ['NCBI', 'ENA661K']

SCHEMA_ANNOTATION_RUNS_CHOICES = ['uniprot-proteomes', 'genbank', 'match-schemas']

# Default values for Schema Refinery's DownloadAssemblies module
FILTERING_CRITERIA = ['abundance', 'genome_size', 'size_threshold',
                      'max_contig_number', 'known_st', 'any_quality',
                      'ST_list_path', 'assembly_level', 'reference',
                      'assembly_source', 'file_to_include', 'verify_status',
                      'exclude_atypical']

ASSEMBLY_LEVELS = ['chromosome', 'complete', 'contig', 'scaffold']

ASSEMBLY_SOURCES = ['RefSeq', 'GenBank', 'all']

FILE_EXTENSIONS = ['genome', 'rna', 'protein',
                   'cds', 'gff3', 'gtf',
                   'gbff', 'seq-report', 'none']

CRITERIA_ERRORS = {'taxon': 
                        ['taxon: must be a valid taxon ID.',
                        [str, None, None, None, None]],
                    'input_table': 
                        ['input_table: must be a valid path to a file.',
                        [str, None, None, True, None]],
                    'abundance':
                        ['abundance: must be float between 0.0 and 1.0.',
                        [float, 0.0, 1.0, None, None]],
                   'genome_size':
                        ['genome_size: must be int greater than 0.',
                        [int, 0, None, None, None]],
                   'size_threshold':
                        ['size_threshold: must be float between 0.0 and 1.0.',
                        [float, 0.0, 1.0, None, None]],
                   'max_contig_number':
                        ['max_contig_number: must be int greater than 0.',
                        [int, 1, None, None, None]],
                   'known_st':
                        ['known_st: must be True or False.',
                        [bool, None, None, None, None]],
                   'any_quality':
                        ['any_quality: must be True or False.',
                        [bool, None, None, None, None]],
                   'ST_list_path':
                        ['ST_list_path: must be a valid path to a file.',
                        [str, None, None, True, None]],
                   'assembly_level':
                        ['assembly_level: one or more of the following separated by a comma: chromosome,complete,contig,scaffold',
                        [str, None, None, None, ASSEMBLY_LEVELS]],
                   'reference':
                        ['reference: must be True or False.',
                        [bool, None, None, None, None]],
                   'assembly_source':
                        ['assembly_source: one of the following: RefSeq,GenBank,all',
                        [str, None, None, None, ASSEMBLY_SOURCES]],
                   'file_to_include':
                        ['file_to_include: one or more of the following separated by comma: genome,rna,protein,cds,gff3,gtf,gbff,seq-report,none',
                        [str, None, None, None, FILE_EXTENSIONS]],
                   'verify_status':
                        ['verify_status: bool or None',
                        [bool, None, None, None, None]],
                   'exclude_atypical':
                        ['exclude_atypical: bool or None',
                        [bool, None, None, None, None]]}

EBI_FTP = 'http://ftp.ebi.ac.uk'
# FTP paths and local location for files needed to download assemblies
ASSEMBLY_FTP_PATH = 'http://ftp.ebi.ac.uk/pub/databases/ENA2018-bacteria-661k/sampleid_assembly_paths.txt'
ASSEMBLY_METADATA_PATH = 'https://figshare.com/ndownloader/files/26578601'
FTP_HASH_FILE = 'http://ftp.ebi.ac.uk/pub/databases/ENA2018-bacteria-661k/checklist.chk'

# IdentifyingSpuriousLoci module

CLASSES_OUTCOMES: Tuple[str, ...] = ('1a', '1b', '2a', '3a', '2b', '1c', '3b', '4a', '4b', '4c', '5')
