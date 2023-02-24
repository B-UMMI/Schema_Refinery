#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
"""


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

CRITERIA_ERRORS = {'abundance':
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
