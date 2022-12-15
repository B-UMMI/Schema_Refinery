#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Purpose
-------

This script checks if the set of alleles for each gene/locus file
given as input respects the specified BSR and length thresholds
applied by chewBBACA's algorithm during schema creation or allele
calling.

Considering that the algorithm is using a BSR >= 0.6...
Every allele that is added to a gene/locus file should have a
BSR >= 0.6 with at least one of the representative alleles of that
gene/locus and should also have a length value that is within the
sequence length mode interval calculated for the  gene/locus. Alleles
from the same gene/locus do not have to always share a BSR >= 0.6  or
be within the sequence length mode interval to be considered valid.
Alleles from a gene/locus do not have to share a BSR >= 0.6 with every
representative allele of that gene/locus. Some alleles might be an ASM
or ALM when compared with other alleles from the same gene/locus but
they will not be an ASM or ALM compared to all representative alleles.
Every allele must be in the accepted sequence length mode range when
compared with at least a representative allele.

For each gene/locus, the full set of alleles is used to create a BLASTp
database. The protein sequences of the representative alleles are BLASTed
against that BLASTp database. BLAST outputs results in tabular format and
the information for each match is completed with additional information,
such as the BLAST Score Ratio and the self raw score from BLAST.

A file with the list of all genes/loci, the number of representative
alleles per gene/locus, the number of alleles that had a valid match with
at least one representative and the alleles that did not have a valid match
is created, as well as a file with the list of all matches between alleles
of the same gene/locus that had a BSR < 0.6.

This script requires that 'Biopython' be installed within the Python 3
environment that is used to run this script.

The main function of this module can be imported and used to obtain the
same results as when the module is called from the CMD.
"""


import os
import sys
import time
import argparse

from Bio import SeqIO

import schema_validation_functions as svf


def main(schema_directory, output_directory, blast_score_ratio, blast_threads):
    """ Determines low-BSR matches between alleles of the same locus
        according to the BLAST Score Ratio value given as threshold.

        Parameters
        ----------
        schema_directory : str
            Path to the schema's directory.
        output_directory : str
            Path to the directory that will store output files.
        blast_score_ratio : float
            BLAST Score Value threshold that will be used to determine
            if the alleles from each locus are represented by at least
            one representative sequence and the level of similarity
            between alleles of the same locus.
        blast_threads : str
            Number of threads for BLASTp.

        Returns
        -------
        Writes a file with matches between alleles of the same locus
        that have a BLAST Score Ratio lower than the defined value.
        Writes another file that shows if each locus passed based on
        the criterion that each allele in a locus file has to have at
        least one match with BSR greater than the defined value with
        one of the representative alleles of that locus.
    """

    start = time.time()
    print('\nImporting and processing schemas...')

    # get list of FASTA files in main schema directory
    loci_files = {file.split('.fasta')[0]: os.path.join(schema_directory, file)
                  for file in os.listdir(schema_directory)
                  if '.fasta' in file}

    # get list of FASTA files in 'short' directory (representative sequences)
    short_directory = os.path.join(schema_directory, 'short')
    # discard files with values of self BSR
    representative_files = {file.split('_short.fasta')[0]: os.path.join(short_directory, file)
                            for file in os.listdir(short_directory)
                            if 'bsr' not in file}

    # create working directory
    if not os.path.exists(output_directory):
        os.mkdir(output_directory)
        print('Created {0} directory to store output '
              'files.\n'.format(output_directory))
    else:
        print('{0} already exists, moving on...\n'.format(output_directory))

    # create dictionary with the DNA sequences of all loci alleles
    fasta_directory = os.path.join(output_directory, 'fastas')
    os.mkdir(fasta_directory)
    reassigned_loci = svf.schema_dictionary(loci_files,
                                            fasta_directory)
    translated_loci = svf.translate_schema(reassigned_loci,
                                           fasta_directory)

    protein_concat = os.path.join(fasta_directory, 'protein_concat.fasta')
    protein_concat = svf.concatenate_files(list(translated_loci.values()),
                                           protein_concat)

    # same steps for the representative sequences
    reassigned_reps = svf.schema_dictionary(representative_files,
                                            fasta_directory)
    translated_representatives = svf.translate_schema(reassigned_reps,
                                                      fasta_directory)

    # create files with recids to BLAST
    ids_files = {}
    for k, v in translated_loci.items():
        recids = [rec.id for rec in SeqIO.parse(v, 'fasta')]
        ids_file = os.path.join(fasta_directory, k+'_ids.txt')
        with open(ids_file, 'w') as outfile:
            outfile.write('\n'.join(recids)+'\n')
        ids_files[k] = ids_file

    # create BLASTdb
    blastdb = os.path.join(fasta_directory, 'protein_concat_blastdb')
    stderr = svf.make_blast_db('makeblastdb', protein_concat, blastdb, 'prot')

    # create list of inputs to distribute with multiprocessing
    # for BLASTp
    blast_files = {k: os.path.join(output_directory, k+'_blastout.tsv')
                   for k in translated_representatives}
    blast_inputs = [['blastp', blastdb, v,
                     blast_files[k],
                     1, 1, ids_files[k],
                     svf.run_blast]
                    for k, v in translated_representatives.items()]

    # BLAST all against all for each locus
    print('Starting BLASTp...')

    blast_stderr = svf.map_async_parallelizer(blast_inputs,
                                              svf.function_helper,
                                              blast_threads,
                                              show_progress=True)

    blast_stderr = svf.flatten_list(blast_stderr)
    if len(blast_stderr) > 0:
        sys.exit(blast_stderr)

    print('Finished BLASTp...')

    # read BLAST results
    blast_results_per_locus = {}
    for locus, file in blast_files.items():
        blast_results_per_locus[locus] = svf.read_blast_tabular(file)
        #os.remove(file)

    # get represetatives
    rep_results = {l[0]: l[2]
                   for k, v in blast_results_per_locus.items()
                   for l in v if l[0] == l[1]}

    alleles_results = {k: [l for l in v if l[0] != l[1]]
                       for k, v in blast_results_per_locus.items()}

    # compute BSR
    for k, v in alleles_results.items():
        for r in v:
            r.append(float(r[2])/float(rep_results[r[0]]))

    # get alleles without BSR > 0.6 with any representative
    low_bsr = {}
    for k, v in alleles_results.items():
        current = {}
        for l in v:
            current.setdefault(l[1], []).append(l[-1])
        low = [(k, v)
               for k, v in current.items()
               if all([r < 0.6 for r in v]) is True]
        if len(low) > 0:
            low_bsr[k] = low

    # write results
    low_bsr_lines = [[(r[0]) for r in v] for k, v in low_bsr.items()]
    low_bsr_lines = svf.flatten_list(low_bsr_lines)
    low_bsr_lines = ['{0}\t{1}'.format(a, str(a in rep_results)) for a in low_bsr_lines]

    low_bsr_file = os.path.join(output_directory, 'low_bsr.tsv')
    with open(low_bsr_file, 'w') as outfile:
        outfile.write('\n'.join(low_bsr_lines)+'\n')

    end = time.time()
    delta = end - start
    print('Done!\nTook {0:.2g} seconds.\n'.format(delta))


def parse_arguments():

    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument('-s', '--schema-directory', type=str,
                        required=True, dest='schema_directory',
                        help='Path to the schema\'s directory or to a file '
                             'that has full paths for the files that should '
                             'be used, one per line.')

    parser.add_argument('-o', '--output-directory', type=str,
                        required=True, dest='output_directory',
                        help='The directory where the output files will be '
                             'saved (will create the directory if it does not '
                             'exist).')

    parser.add_argument('-bsr', '--blast-score-ratio', type=float,
                        required=False, dest='blast_score_ratio',
                        default=0.6,
                        help='The BLAST Score Ratio value that '
                             'will be used to determine if the '
                             'schema is valid (default=0.6).')

    parser.add_argument('-bt', '--blast-threads', type=int,
                        required=False, dest='blast_threads',
                        default=1,
                        help='The number of threads to pass as argument '
                             'to BLASTp (default=1).')

    args = parser.parse_args()

    return args


if __name__ == "__main__":

    args = parse_arguments()
    main(**vars(args))
