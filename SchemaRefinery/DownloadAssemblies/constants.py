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
                   'gbff', 'seq-report']

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
                       ['known_st: must bool or None (Empty or None is False).',
                        [bool, None, None, None, None]],
                   'any_quality':
                       ['any_quality: must bool or None (Empty or None is True).',
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
                       ['verify_status: bool or None (Empty or None is False)',
                        [bool, None, None, None, None]],
                   'exclude_atypical':
                       ['exclude_atypical: bool or None (Empty or None is True)',
                        [bool, None, None, None, None]]}

EBI_FTP = 'http://ftp.ebi.ac.uk'
# FTP paths and local location for files needed to download assemblies
ASSEMBLY_FTP_PATH = 'http://ftp.ebi.ac.uk/pub/databases/ENA2018-bacteria-661k/sampleid_assembly_paths.txt'
ASSEMBLY_METADATA_PATH = 'https://figshare.com/ndownloader/files/26578601'
FTP_HASH_FILE = 'http://ftp.ebi.ac.uk/pub/databases/ENA2018-bacteria-661k/checklist.chk'


"""
Modify_schema
"""

ACCEPTED_COMMANDS = ['merge', 'split', 'remove']

INPUT_ERRORS = {'merge':
                    ['merge: must contain string values.',
                    [str, None, None, None, None]],
                'split_id':
                    ['split: must contain a string with the new id.',
                    [str, None, None, None, None]],
                'split_size':
                    ['split: must have the sizes superior to 0.',
                    [int, 0, None, None, None]],
                'remove':
                    ['split: must contain a string.',
                    [str, None, None, None, None]]}

SPLIT_MUST_BE_ODD = "split: must be with the following format: split    loci5   new_loci_x   100-300   new_loci_y   301-500."

SPLIT_MISSING_MINUS = "split: must have '-' between two values."

SPLIT_VALUE_MUST_BE_INT = "split: must have integer number separated by '-'"


"""
ncbi_linked_ids
"""

# regex expressions to identify identifier type
DATABASE_PATTERNS = {'biosample': 'SAM[E|D|N][A-Z]?[0-9]+',
                     'bioproject': 'PRJ[E|D|N][A-Z][0-9]+',
                     'sra': '[E|D|S]RR[0-9]{6,}',
                     'refseq': 'GCF_[0-9]{9}.[0-9]+',
                     'genbank': 'GCA_[0-9]{9}.[0-9]+'}