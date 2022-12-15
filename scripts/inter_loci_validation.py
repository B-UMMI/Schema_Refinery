#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Purpose
-------

This script verifies if a set of FASTA files corresponding to
different loci respect the specified BSR and length thresholds
applied by chewBBACA's algorithm during schema creation or
allele calling.

Considering that the algorithm is using a BSR >= 0.6...
Representative alleles from different loci should not match with
a BSR >= 0.6. Such cases are considered to be paralogs. A schema
created by the CreateSchema process should contain only a
representative allele for each locus/gene in the schema and each
representative allele must not share a BSR >= 0.6 with any of the
representative alleles for the other loci/genes. The AlleleCall
process adds new alleles to a gene/locus if those new alleles share
a BSR >= 0.6 with any representative allele for that locus/gene and
are contained in the allele sequence lenght mode interval. However,
the new allele does not have to share a BSR >= 0.6 with all alleles
for that locus/gene and may also share a BSR >= 0.6 with
non-representative alleles of other loci/genes (and with representatives
if it is non-representative).

This script outputs files with relevant information that help identify
problematic cases in which representative alleles from different loci/genes
match with a BSR >= 0.6. This can be helpful in verifying if a new schema
created with the CreateSchema process has representative alleles that are
all BSR < 0.6. If not, some alleles that should be considered as belonging
to the same locus/gene were separated into different loci/genes and will be
listed as paralogs in the AlleleCall process, possibly leading to the removal
and loss of those genes that might have stayed in the schema otherwise. For
schemas that have been populated with additional alleles through the
AlleleCall process it is relevant to verify if all representative alleles
for different loci still share a BSR < 0.6 and to list all matches between
alleles of different loci that have a BSR >= 0.6 (representative or not).
Information about sequence length is also used to determine if matches with
BSR >= 0.6 correspond to ASMs or ALMs.

The full set of alleles from all genes/loci is used to create a BLASTp
database. The same protein sequences used to create the BLASTp database are
BLASTed against the database to find all matches between alleles of different
loci. BLAST outputs results in tabular format and the information for each
match is completed with additional information, as the BLAST Score Ratio,
the self raw score from BLAST and if any of the matched sequences is a
representative.

A file with the list of all gene/loci and information about if they are
valid is created, as well as a file that contains the information about
the matches between alleles of different genes/loci that have a BSR >= 0.6.

This script requires that 'Biopython' be installed within the Python 3
environment that is used to run this script.

The main function of this module can be imported and used to obtain the
same results as when the module is called from the CMD.

Code documentation
------------------
"""


import os
import csv
import sys
import time
import argparse

from Bio import SeqIO

# these modules must be in the same directory
import schema_validation_functions as svf


def exportListToTSVFile(filename: str, lines_list: list):

    with open(filename, 'w') as csvfile:
        # creating a csv writer object
        csvwriter = csv.writer(csvfile, delimiter='\t')

        # writing the data rows
        csvwriter.writerows(lines_list)


def main(schema_directory, output_directory, blast_score_ratio, blast_threads):
    """Identify alleles from different loci that are highly similar.

    Parameters
    ----------
    schema_directory : str
        path to the schema with all loci files.
    output_directory : str
        path/name of folder that will store output files.
    blast_score_ratio : float
        threshold BLAST Score Ratio value that will be used
        to determine if alleles from different loci share
        greater than expected similarity.
    num_threads : str
        number of threads for BLASTp searches.

    Returns
    -------
    Writes a file with matches between alleles of different
    loci that have a BLAST Score Ratio greater than the defined
    value. Writes another file that shows, for each locus, if
    there were high-BSR matches with alleles from other loci,
    with which alleles and if the loci passed based on the
    criterion that representative alleles from different loci
    should not have matches with a BSR greater than the defined
    value.
    """
    start = time.time()
    print('\nImporting and processing loci...')

    # get list of FASTA files in main schema directory
    loci_files = {file.split('.fasta')[0]: os.path.join(schema_directory, file)
                  for file in os.listdir(schema_directory)
                  if '.fasta' in file}

    # get list of FASTA files in 'short' directory (representative sequences)
    short_directory = os.path.join(schema_directory, 'short')
    # discard files with values of self BSR
    representative_files = {file.split('_short.fasta')[0]: os.path.join(short_directory, file)
                            for file in os.listdir(short_directory)
                            if file.endswith('fasta')}

    # create working directory
    if not os.path.exists(output_directory):
        os.mkdir(output_directory)
        print('Created {0} directory to store output '
              'files.\n'.format(output_directory))
    else:
        print('{0} already exists, moving on...\n'.format(output_directory))

    # create files with translated alleles for all loci
    fasta_directory = os.path.join(output_directory, 'fastas')
    os.mkdir(fasta_directory)
    reassigned_loci = svf.schema_dictionary(loci_files,
                                            fasta_directory)
    translated_loci = svf.translate_schema(reassigned_loci,
                                           fasta_directory)

    # same steps for the representative sequences
    reassigned_reps = svf.schema_dictionary(representative_files,
                                            fasta_directory)
    translated_representatives = svf.translate_schema(reassigned_reps,
                                                      fasta_directory)

    # concatenate all representatives
    rep_concat = os.path.join(fasta_directory, 'rep_concat.fasta')
    rep_concat = svf.concatenate_files(list(translated_representatives.values()),
                                       rep_concat)

    # BLAST representatives to determine self-score
    # It is important to BLAST representatives separately because the
    # maximum number of reported alignments might be reached if we align
    # against the full schema and we might not get the self-alignment for
    # all representatives
    rep_blastdb = os.path.join(fasta_directory, 'rep_concat_blastdb')
    stderr = svf.make_blast_db('makeblastdb', rep_concat, rep_blastdb, 'prot')

    rep_blast_inputs = []
    for k, v in translated_representatives.items():
        locus_id = k
        locus_blastout = os.path.join(output_directory, locus_id+'_blastout.tsv')
        # create file with representative identifiers
        locus_ids_file = os.path.join(output_directory, locus_id+'_ids.txt')
        locus_ids = [rec.id for rec in SeqIO.parse(v, 'fasta')]
        with open(locus_ids_file, 'w') as outfile:
            outfile.write('\n'.join(locus_ids)+'\n')

        rep_blast_inputs.append(['blastp', rep_blastdb, v, locus_blastout, 1, 1, locus_ids_file, svf.run_blast])

    # BLAST all against all for each locus
    print('BLASTing representative sequences to determine self-score...')
    blast_stderr = svf.map_async_parallelizer(rep_blast_inputs,
                                              svf.function_helper,
                                              blast_threads,
                                              show_progress=True)

    blast_stderr = svf.flatten_list(blast_stderr)
    if len(blast_stderr) > 0:
        sys.exit(blast_stderr)

    # concatenate all files with BLAST results
    rep_blastout_files = [os.path.join(output_directory, file)
                          for file in os.listdir(output_directory)
                          if file.endswith('blastout.tsv') is True]

    # concatenate all BLAST output files
    rep_blast_output = os.path.join(output_directory, 'reps_blastout_concat.tsv')
    rep_blast_output = svf.concatenate_files(rep_blastout_files, rep_blast_output)

    # read BLAST results
    representative_results = svf.read_blast_tabular(rep_blast_output)

    # get representatives self-score
    rep_results = {l[0]: l[2] for l in representative_results if l[0] == l[1]}

    # create the FASTA file with all protein sequences for the BLASTp database
    protein_concat = os.path.join(fasta_directory, 'protein_concat.fasta')
    protein_concat = svf.concatenate_files(list(translated_loci.values()),
                                           protein_concat)

    # create BLASTdb
    blastdb = os.path.join(fasta_directory, 'protein_concat_blastdb')
    stderr = svf.make_blast_db('makeblastdb', protein_concat, blastdb, 'prot')

    # create list of inputs to distribute with multiprocessing
    # for BLASTp
    blast_files = {k: os.path.join(output_directory, k+'_blastout_main.tsv')
        for k in translated_representatives}
    blast_inputs = []
    for i in rep_blast_inputs:
        locus_blastout = i[3].replace('.tsv', '_main.tsv')
        blast_inputs.append(['blastp', blastdb, i[2], locus_blastout, svf.run_blast])

    # BLAST all against all for each locus
    print('BLASTing representative sequences against schema...')
    blast_stderr = svf.map_async_parallelizer(blast_inputs,
                                              svf.function_helper,
                                              blast_threads,
                                              show_progress=True)

    blast_stderr = svf.flatten_list(blast_stderr)
    if len(blast_stderr) > 0:
        sys.exit(blast_stderr)

    # read BLAST results
    blast_results_per_locus = {}
    for locus, file in blast_files.items():
        blast_results_per_locus[locus] = svf.read_blast_tabular(file)

    alleles_results = {k: [l for l in v if l[0] != l[1]]
                       for k, v in blast_results_per_locus.items()}

    # compute BSR
    for k, v in alleles_results.items():
        for r in v:
            r.append(str(float(r[2])/float(rep_results[r[0]])))

    high_bsr = []
    for k, v in alleles_results.items():
        current = [r
                   for r in v
                   if '_'.join(r[0].split('_')[0:-1]) != '_'.join(r[1].split('_')[0:-1])
                   and float(r[-1]) >= 0.6]
        if len(current) > 0:
            high_bsr.extend(current)

    # index Fasta
    index = SeqIO.index(protein_concat, 'fasta')
    for r in high_bsr:
        query = index.get(r[0])
        query = (query.id, str(query.seq))
        subject = index.get(r[1])
        subject = (subject.id, str(subject.seq))
        # determine if ASM or ALM
        default = 'IN_MODE'
        if len(query[1]) < (len(subject[1])-(0.2*len(subject[1]))):
            default = 'ASM'
        elif len(query[1]) > (len(subject[1])+(0.2*len(subject[1]))):
            default = 'ALM'

        r.append(default)

    # save results
    output_file = os.path.join(output_directory, 'high_bsr.tsv')
    exportListToTSVFile(output_file, high_bsr)

    end = time.time()
    delta = end - start
    print('Done!\nTook {0} seconds.\n'.format(delta))
    print('Results available at {0}'.format(output_file))


def parse_arguments():

    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument('-s', '--schema-directory', type=str,
                        required=True, dest='schema_directory',
                        help='Path to the schema\'s directory.')

    parser.add_argument('-o', '--output-directory', type=str,
                        required=True, dest='output_directory',
                        help='The directory where the output files will be '
                             'saved (will create the directory if it does not '
                             'exist).')

    parser.add_argument('-bsr', '--blast-score-ratio', type=float,
                        required=False, dest='blast_score_ratio',
                        default=0.6,
                        help='The BLAST Score Ratio value that '
                             'will be used to compare '
                             'loci and identify paralogs.')

    parser.add_argument('-bt', '--blast-threads', type=int,
                        required=False, dest='blast_threads',
                        default=1,
                        help='The number of threads to pass '
                             'to BLASTp.')

    args = parser.parse_args()

    return args


if __name__ == "__main__":

    args = parse_arguments()
    main(**vars(args))
