#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
AUTHOR

    Rafael Mamede
    github: @rfm-targa

DESCRIPTION

    This script checks if the set of alleles for each gene/locus file given as 
    input respects the specified BSR and length thresholds applied by chewBBACA's
    algorithm during schema creation or allele calling.

    Considering that the algorithm is using a BSR >= 0.6...
    Every allele that is added to a gene/locus file should have a BSR >= 0.6 with at 
    least one of the representative alleles of that gene/locus and should also have a
    length value that is within the sequence length mode interval calculated for the 
    gene/locus. Alleles from the same gene/locus do not have to always share a BSR >= 0.6 
    or be within the sequence length mode interval to be considered valid. Alleles from a 
    gene/locus do not have to share a BSR >= 0.6 with every representative allele of that 
    gene/locus. Some alleles might be an ASM or ALM when compared with other alleles from 
    the same gene/locus but they will not be an ASM or ALM compared to all representative 
    alleles. Every allele must be in the accepted sequence length mode range when compared 
    with at least a representative allele.

    For each gene/locus, the full set of alleles is used to create a BLASTp database. 
    The protein sequences of the representative alleles are BLASTed against that 
    BLASTp database. BLAST outputs results in tabular format and the information 
    for each match is completed with additional information, such as the BLAST Score 
    Ratio and the self raw score from BLAST.

    A file with the list of all genes/loci, the number of representative alleles per 
    gene/locus, the number of alleles that had a valid match with at least one 
    representative and the alleles that did not have a valid match is created, as well 
    as a file with the list of all matches between alleles of the same gene/locus that 
    had a BSR < 0.6.
    
    This script requires that 'Biopython' be installed within the Python 3 environment 
    that is used to run this script.
    
    The main function of this module can be imported and used to obtain the same 
    results as when the module is called from the CMD.
"""


import os
import sys
import time
import argparse
from multiprocessing import Pool

import schema_validation_functions as svf


def main(schema_full, schema_short, output_directory, blast_score_ratio, num_threads, files_prefix):
    """ Determines low-BSR matches between alleles of the same locus according to the 
        BLAST Score Ratio value given as threshold.

        Args:
            schema_full (str): path to the schema with all loci files.
            schema_short (str): path to the short schema with all representative
            alleles files.
            blast_score_ratio (float): BLAST Score Value threshold that will be used
            to determine if the alleles from each locus are represented by at least
            one representative sequence and the level of similarity between alleles of
            the same locus.
            num_threads (str): number of threads for BLASTp searches.
            output_directory (str): path/name of folder that will store output files.
            files_prefix (str): prefix for the output files.

        Returns:
            Writes a file with matches between alleles of the same locus that have
            a BLAST Score Ratio lower than the defined value. Writes another file
            that shows if each locus passed based on the criterion that each allele
            in a locus file has to have at least one match with BSR greater than the
            defined value with one of the representative alleles of that locus.
    """

    start = time.time()
    print('\nImporting and processing schemas...')

    # get list of FASTA files in main schema directory
    full_base_schema, schema_directory = svf.identify_input_type(schema_full)
    # remove 'short' and 'protogenome' directories from list of files
    full_base_schema = [file for file in full_base_schema if file not in ['short', 'protogenome']]

    # get list of FASTA files in 'short' directory (representative sequences)
    short_base_schema, short_directory = svf.identify_input_type(schema_short)
    # discard files with values of self BSR
    short_base_schema = [file for file in short_base_schema if 'bsr' not in file]
    print(short_base_schema)
    # create working directory
    if not os.path.exists(output_directory):
        os.mkdir(output_directory)
        print('Created {0} directory to store output files.\n'.format(output_directory))
    else:
        print('{0} already exists, moving on...\n'.format(output_directory))

    # create dictionary with the DNA sequences of all loci alleles
    full_dna_schema, genes_ids, alleles_ids = svf.schema_dictionary(full_base_schema, schema_directory)

    # create a new dictionary with the protein sequences
    full_protein_schema = svf.translate_schema(full_dna_schema)

    # same steps for the representative sequences
    print(short_directory)
    short_dna_schema = svf.short_schema_dict(short_base_schema, short_directory,
                                             genes_ids, alleles_ids)
    short_protein_schema = svf.translate_schema(short_dna_schema)

    representatives_identifiers = svf.determine_representatives(short_protein_schema)
    sys.exit(0)
    # create FASTA lines with all proteins for each gene and save in files
    full_protein_lines = svf.schema_fasta_lines(full_protein_schema)

    full_fasta_paths = svf.write_per_locus_fasta_file(full_protein_lines, output_directory, '')

    # create FASTA lines with the representative protein sequences and save in files
    short_protein_lines = svf.schema_fasta_lines(short_protein_schema)
    short_fasta_paths = svf.write_per_locus_fasta_file(short_protein_lines, output_directory, '_short')

    # create list of inputs to distribute with multiprocessing
    # for BLASTp
    blast_inputs = []
    for p in range(len(full_fasta_paths)):
        blast_inputs.append([full_fasta_paths[p], short_fasta_paths[p], output_directory, 1])

    # BLAST all against all for each locus
    print('Starting BLASTp...')

    blast_files = []
    pool = Pool(processes=num_threads)
    r = pool.map_async(svf.within_loci_blast, blast_inputs, callback=blast_files.append)
    r.wait()

    print('Finished BLASTp...')

    # read BLAST results
    blast_results_per_locus = []
    for file in blast_files[0]:
        blast_results_per_locus.append(svf.read_blast_tabular(file))
        os.remove(file)

    within_loci_results = []
    for locus_matches in blast_results_per_locus:
        blast_results_rows = svf.create_blast_results_dictionary(locus_matches)

        within_loci_results.append(blast_results_rows)

    print('Determining BLASTp matches information and matches with BLAST Score Ratio < {0}...\n'.format(str(blast_score_ratio)))

    # determine all info for all matches
    informative_matches = []
    bsr_distances = []
    for locus in within_loci_results:
        try:
            complete_info_match, distance_matrix = svf.complete_loci_matches_info(locus, full_protein_schema,
                                                                                  representatives_identifiers)
        except:
            print(locus, representatives_identifiers)
            sys.exit(1)
        if len(complete_info_match) > 0:
            informative_matches.append(complete_info_match)
            bsr_distances.append(distance_matrix)

    files_header = ('Query,Subject,Identity,Alignment_length,Query_start,Query_end,'
                    'Subject_start,Subject_end,Bitscore,Raw_score,Self_score,BSR,'
                    'Query_length,Subject_length,Locus_length_mode,Match_type,'
                    'Representative\n\n')

    all_matches_lines, low_bsr_lines = svf.create_within_loci_lines(informative_matches, blast_score_ratio,
                                                                    files_header, genes_ids, alleles_ids)

    # save files with low-BSR matches between alleles of same locus
    low_bsr_matches_outfile = '{0}/{1}_low_bsr_matches.csv'.format(os.path.abspath(output_directory), files_prefix)
    svf.write_list_to_file(low_bsr_lines, low_bsr_matches_outfile)

    # determine alleles that passed the criteria (all alleles must share one match with a representative with a BSR > 0.6)
    # get number and set of alleles per locus
    alleles_info = svf.allele_set(full_protein_schema)

    matched_representative, not_matched = svf.validate_same_locus(informative_matches, alleles_info,
                                                                  blast_score_ratio, 11, genes_ids)

    matched_header = 'Locus,Number_Representatives,Number_alleles,Matched_alleles,Validate,Unmatched_Alleles\n'

    if len(matched_representative) > 0:
        matched_file = '{0}/{1}_intra_loci.csv'.format(os.path.abspath(output_directory), files_prefix)
        svf.write_matches_validation(matched_representative, matched_file, matched_header)

    num_of_matches = len(informative_matches)

    print(('Found {0}/{1} cases in which alleles from the same locus do not have a BSR >= {2} with at least one '
           'representative sequence of the locus.\n').format(not_matched, num_of_matches, blast_score_ratio))

    num_of_low_bsr = len([line for line in low_bsr_lines[1:] if len(line) > 1])
    print('Found {0} cases in which the alignment between a representative sequence and an allele has a '
          'BSR < {1}.'.format(num_of_low_bsr, blast_score_ratio))

    number_of_rep = len([1 for line in low_bsr_lines if ',11\n' in line])
    number_of_asm = len([1 for line in low_bsr_lines if 'ASM' in line])
    number_of_alm = len([1 for line in low_bsr_lines if 'ALM' in line])

    print('\t{0} with both sequences as representatives;\n\t{1} are ASM cases;\n\t{2} are ALM cases.\n'.format(number_of_rep, number_of_asm, number_of_alm))

    print('Wrote evaluation results to: \n{0} \n{1}\n'.format(low_bsr_matches_outfile, matched_file))

    end = time.time()
    delta = end - start
    print('Done!\nTook {0:.2g} seconds.\n'.format(delta))


def parse_arguments():

    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument('-fs', '--full_schema', type=str, required=True, dest='full_schema_files',
                        help='Path to the directory with the schema FASTA files that have all alleles '
                             'or to a file that has full paths for the files that should be used, '
                             'one per line.')

    parser.add_argument('-ss', '--short_schema', type=str, required=True, dest='short_schema_files',
                        help='Path to the directory with the schema FASTA files that only have the '
                             'representative alleles or to a file that has full paths for the files '
                             'that should be used, one per line.')

    parser.add_argument('-o', '--output_directory', type=str, required=True, dest='output_directory',
                        help='The directory where the output files will be saved (will create the directory if '
                             'it does not exist).')

    parser.add_argument('-bsr', '--blast_score_ratio', type=float, required=False, dest='blast_score_ratio',
                        default=0.6, help='The BLAST Score Ratio value that will be used to determine if the '
                                          'schema is valid (default=0.6).')

    parser.add_argument('-bt', '--blast_threads', type=int, required=False, dest='blast_threads',
                        default=1, help='The number of threads to pass as argument to BLASTp (default=1).')

    parser.add_argument('-fp', '--files_prefix', type=str, required=False, dest='output_files_prefix',
                        default='ilv', help='A prefix for the output files names (default=ilv).')

    args = parser.parse_args()

    return [args.full_schema_files, args.short_schema_files, args.output_directory, 
            args.blast_score_ratio, args.blast_threads, args.output_files_prefix]


if __name__ == "__main__":

    args = parse_arguments()
    main(args[0], args[1], args[2], args[3], args[4], args[5])
