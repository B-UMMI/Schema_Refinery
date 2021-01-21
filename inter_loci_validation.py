#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
AUTHOR

    Rafael Mamede
    github: @rfm-targa

DESCRIPTION

    This script verifies if a set of FASTA files corresponding to different loci 
    respect the specified BSR and length thresholds applied by chewBBACA's algorithm 
    during schema creation or allele calling.

    Considering that the algorithm is using a BSR >= 0.6...
    Representative alleles from different loci should not match with a BSR >= 0.6. 
    Such cases are considered to be paralogs. A schema created by the CreateSchema 
    process should contain only a representative allele for each locus/gene in the 
    schema and each representative allele must not share a BSR >= 0.6 with any of 
    the representative alleles for the other loci/genes. The AlleleCall process adds 
    new alleles to a gene/locus if those new alleles share a BSR >= 0.6 with any 
    representative allele for that locus/gene and are contained in the allele sequence 
    lenght mode interval. However, the new allele does not have to share a BSR >= 0.6 
    with all alleles for that locus/gene and may also share a BSR >= 0.6 with 
    non-representative alleles of other loci/genes (and with representatives if it is 
    non-representative).

    This script outputs files with relevant information that help identify problematic 
    cases in which representative alleles from different loci/genes match with a 
    BSR >= 0.6. This can be helpful in verifying if a new schema created with the 
    CreateSchema process has representative alleles that are all BSR < 0.6. If not, 
    some alleles that should be considered as belonging to the same locus/gene were 
    separated into different loci/genes and will be listed as paralogs in the AlleleCall 
    process, possibly leading to the removal and loss of those genes that might have 
    stayed in the schema otherwise. For schemas that have been populated with additional 
    alleles through the AlleleCall process it is relevant to verify if all representative 
    alleles for different loci still share a BSR < 0.6 and to list all matches between 
    alleles of different loci that have a BSR >= 0.6 (representative or not). Information 
    about sequence length is also used to determine if matches with BSR >= 0.6 correspond 
    to ASMs or ALMs.

    The full set of alleles from all genes/loci is used to create a BLASTp database. 
    The same protein sequences used to create the BLASTp database are BLASTed against 
    the database to find all matches between alleles of different loci. BLAST outputs 
    results in tabular format and the information for each match is completed with 
    additional information, as the BLAST Score Ratio, the self raw score from BLAST 
    and if any of the matched sequences is a representative.

    A file with the list of all gene/loci and information about if they are valid is 
    created, as well as a file that contains the information about the matches between 
    alleles of different genes/loci that have a BSR >= 0.6.

    This script requires that 'Biopython' be installed within the Python 3 environment 
    that is used to run this script.

    The main function of this module can be imported and used to obtain the same 
    results as when the module is called from the CMD.
"""


import os
import time
import argparse

import schema_validation_functions as svf


def main(schema_full, schema_short, output_directory, blast_score_ratio, num_threads, files_prefix):
    """ Determines BLASTp matches between different loci of a schema with a BSR greater than 
        the defined value.
        
        Args:
            schema_full (str): path to the schema with all loci files.
            schema_short (str): path to the short schema with all representative
            alleles files.
            blast_score_ratio (float): threshold BLAST Score Ratio value that will be used 
            to determine if alleles from different loci share greater than expected similarity.
            num_threads (str): number of threads for BLASTp searches.
            output_directory (str): path/name of folder that will store output files.
            files_prefix (str): prefix for the output files.
            
        Returns:
            Writes a file with matches between alleles of different loci that have 
            a BLAST Score Ratio greater than the defined value. Writes another file 
            that shows, for each locus, if there were high-BSR matches with alleles 
            from other loci, with which alleles and if the loci passed based on the
            criterion that representative alleles from different loci shoud not have 
            matches with a BSR greater than the defined value.
    """

    start = time.time()
    print('\nImporting and processing schemas...')

    # get list of files in schema directory
    full_base_schema, schema_directory = svf.identify_input_type(schema_full)
    full_base_schema = [file for file in full_base_schema if file not in ['short', 'protogenome']]

    short_base_schema, short_directory = svf.identify_input_type(schema_short)
    short_base_schema = [file for file in short_base_schema if 'bsr' not in file]

    # create dictionaries with genes/loci identifiers as keys and DNA or protein
    # sequences as values
    full_dna_schema, genes_ids, alleles_ids = svf.schema_dictionary(full_base_schema, schema_directory)
    full_protein_schema = svf.translate_schema(full_dna_schema)

    short_dna_schema = svf.short_schema_dict(short_base_schema, short_directory,
                                             genes_ids, alleles_ids)
    short_protein_schema = svf.translate_schema(short_dna_schema)

    representatives_identifiers = svf.determine_representatives(short_protein_schema)

    # create directory to store results
    if not os.path.exists(output_directory):
        os.mkdir(output_directory)
        print('Created {0} directory to store output files.\n'.format(output_directory))
    else:
        print('{0} already exists, moving on...\n'.format(output_directory))

    # create the FASTA file with all protein sequences for the BLASTp database
    schema_lines = svf.schema_fasta_lines(full_protein_schema)
    fasta_file = '{0}/{1}_alleles.fasta'.format(output_directory,files_prefix)
    svf.write_all_loci_fasta_file(schema_lines, fasta_file)

    # create blastdb
    print('Creating BLASTp database...')
    blastdb_name = svf.create_blastdb(output_directory, fasta_file, True)

    print('Starting BLASTp...')
    # use BLAST to check if any pair of sequences from diffrent loci have BSR >= 0.6
    blast_out_file = '{0}/{1}_blast_out.tsv'.format(output_directory,files_prefix)
    blast_command = 'blastp -db {0} -query {1} -out {2} -outfmt "6 qseqid sseqid pident length qstart qend sstart send evalue bitscore score" -max_hsps 1 -max_target_seqs 100000 -num_threads {3} -evalue 0.0001'.format(blastdb_name, fasta_file, blast_out_file, num_threads)
    os.system(blast_command)
    print('Finished BLASTp!')

    # read BLASTp results
    blast_results = svf.read_blast_tabular(blast_out_file)

    # create dictionary with blast results organized in a hierarchical way
    blast_results_rows = svf.create_blast_results_dictionary(blast_results)

    betwenn_loci_results = []
    for locus in blast_results_rows:
        betwenn_loci_results.append({locus:blast_results_rows[locus]})

    # determine additional information for each match (self BLAST score, BSR, locus length mode...)

    print('Determining BLASTp matches information and matches with BLAST Score Ratio >= {0}...\n'.format(str(blast_score_ratio)))

    informative_matches = []
    bsr_distances = []
    for locus in betwenn_loci_results:
        complete_info_match, distance_matrix = svf.complete_loci_matches_info(locus, full_protein_schema, representatives_identifiers)
        if len(complete_info_match) > 0:
            informative_matches.append(complete_info_match)
            bsr_distances.append(distance_matrix)

    # Determine cases where alleles from one locus match with alleles of other loci
    different_loci_matches = []
    for locus in informative_matches:
        different_loci_matches.append(svf.inter_loci_matches(locus))

    # write all matches with different alleles and just the ones that have BSR >= 0.6
    files_header = ('Query,Subject,Identity,Alignment_length,Query_start,Query_end,'
                   'Subject_start,Subject_end,Bitscore,Raw_score,Self_score,BSR,'
                   'Query_length,Subject_length,Locus_length_mode,Match_type,Representative\n\n')

    all_matches_lines, high_bsr_lines, passed = svf.validate_different_loci(different_loci_matches, representatives_identifiers, 
                                                                            files_header, blast_score_ratio, genes_ids, alleles_ids)

    # save files with matches
    high_bsr_matches_file = '{0}/{1}_high_bsr_matches.csv'.format(os.path.abspath(output_directory), files_prefix)
    with open(high_bsr_matches_file, 'w') as lbl:
        lbl.writelines(high_bsr_lines)

    # Write cases that passed and cases that did not pass
    loci_verdict_file = '{0}/{1}_inter_loci.csv'.format(os.path.abspath(output_directory), files_prefix)
    with open(loci_verdict_file, 'w') as npl:
        passed_lines = [k+','+v+'\n' for k,v in passed.items()]
        passed_lines = ['Gene,Representative_Matches,Validate\n'] + passed_lines
        npl.writelines(passed_lines)

    num_not_passed_loci = len([line[0] for line in passed_lines if 'NOT_PASSED' in line])
    print(('Found {0} loci with representative sequences that match with representative sequences '
           'of other loci with a BSR >= {1}.\n').format(num_not_passed_loci, blast_score_ratio))

    num_of_high_bsr = len([line for line in high_bsr_lines[1:] if len(line) > 1])
    print(('Found {0} cases in which alleles from different loci have a BLASTp match with a BSR >= {1}.').format(num_of_high_bsr, blast_score_ratio))

    number_of_rep = len([1 for line in high_bsr_lines if ',11\n' in line])
    number_of_asm = len([1 for line in high_bsr_lines if 'ASM' in line])
    number_of_alm = len([1 for line in high_bsr_lines if 'ALM' in line])

    print('\t{0} with both sequences as representatives;\n\t{1} are ASM cases;\n\t{2} are ALM cases.\n'.format(number_of_rep, number_of_asm, number_of_alm))

    print('Evaluation results stored in: \n{0} \n{1}\n'.format(high_bsr_matches_file, loci_verdict_file))

    end = time.time()
    delta = end - start
    print('Done!\nTook {0} seconds.\n'.format(delta))


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
