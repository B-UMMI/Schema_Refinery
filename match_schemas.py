#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Purpose
-------
This script aligns representative alleles in a query schema
against all alleles in a subject schema to determine similar
loci in both schemas.

Code documentation
------------------
"""


import os
import sys
import csv
import argparse
import itertools
import subprocess

from Bio import SeqIO
from Bio.Seq import Seq


def reverse_str(string):
    """ Reverse character order in input string.

    Parameters
    ----------
    string : str
        String to be reversed.

    Returns
    -------
    revstr : str
        Reverse of input string.
    """

    revstr = string[::-1]

    return revstr


def reverse_complement(dna_sequence):
    """ Determines the reverse complement of a DNA sequence.

    Parameters
    ----------
    dna_sequence : str
        String representing a DNA sequence.

    Returns
    -------
    reverse_complement : str
        The reverse complement of the input DNA
        sequence.
    """

    base_complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}

    # convert string into list with each character as a separate element
    bases = list(dna_sequence.upper())

    # determine complement strand
    complement_bases = []
    for base in bases:
        complement_bases.append(base_complement.get(base, base))

    complement_strand = ''.join(complement_bases)

    # reverse strand
    reverse_complement = reverse_str(complement_strand)

    return reverse_complement


def translate_sequence(dna_str, table_id):
    """ Translate a DNA sequence using the BioPython package.

    Parameters
    ----------
    dna_str : str
        DNA sequence.
    table_id : int
        Translation table identifier.

    Returns
    -------
    protseq : str
        Protein sequence created by translating
        the input DNA sequence.
    """

    myseq_obj = Seq(dna_str)
    # sequences must be a complete and valid CDS
    protseq = Seq.translate(myseq_obj, table=table_id, cds=True)

    return protseq


def make_blast_db(input_fasta, output_path, db_type):
    """ Creates a BLAST database.

    Parameters
    ----------
    input_fasta : str
        Path to a Fasta file.
    output_path : str
        Path to the output BLAST database.
    db_type : str
        Type of the database, nucleotide (nuc) or
        protein (prot).
    """

    blastdb_cmd = ['makeblastdb', '-in', input_fasta, '-out', output_path,
                   '-parse_seqids', '-dbtype', db_type]

    makedb_cmd = subprocess.Popen(blastdb_cmd,
                                  stdout=subprocess.PIPE,
                                  stderr=subprocess.PIPE)

    stdout, stderr = makedb_cmd.communicate()
    print(stdout, stderr)

    makedb_cmd.wait()


def run_blast(blast_path, blast_db, fasta_file, blast_output,
              max_hsps=1, threads=1, ids_file=None, blast_task=None,
              max_targets=None):
    """ Executes BLAST.

    Parameters
    ----------
    blast_path : str
        Path to the BLAST executable.
    blast_db : str
        Path to the BLAST database.
    fasta_file : str
        Path to the Fasta file that contains the sequences
        to align against the database.
    blast_output : str
        Path to the output file.
    max_hsps : int
        Maximum number of High-Scoring Pairs.
    threads : int
        Number of threads passed to BLAST.
    ids_file : path
        Path to a file with the identifiers of the sequences
        to align against. Used to specify the database sequences
        we want to align against.
    blast_task : str
        BLAST task. Allows to set default parameters for a specific
        type of search.
    max_targets : int
        Maximum number of targets sequences to align against.

    Returns
    -------
    stderr : list
        List with the warnings/errors reported by BLAST.
    """

    blast_args = [blast_path, '-db', blast_db, '-query', fasta_file,
                  '-out', blast_output, '-outfmt', '6 qseqid sseqid score',
                  '-max_hsps', str(max_hsps), '-num_threads', str(threads),
                  '-evalue', '0.001']

    if ids_file is not None:
        blast_args.extend(['-seqidlist', ids_file])
    if blast_task is not None:
        blast_args.extend(['-task', blast_task])
    if max_targets is not None:
        blast_args.extend(['-max_target_seqs', str(max_targets)])

    blast_proc = subprocess.Popen(blast_args,
                                  stdout=subprocess.PIPE,
                                  stderr=subprocess.PIPE)

    stderr = blast_proc.stderr.readlines()

    return stderr


def read_tabular(input_file, delimiter='\t'):
    """ Read a TSV file.

    Parameters
    ----------
    input_file : str
        Path to a tabular file.
    delimiter : str
        Delimiter used to separate file fields.

    Returns
    -------
    lines : list
        A list with a sublist per line in the input file.
        Each sublist has the fields that were separated by
        the specified delimiter.
    """

    with open(input_file, 'r') as infile:
        reader = csv.reader(infile, delimiter=delimiter)
        lines = [line for line in reader]

    return lines


def flatten_list(list_to_flatten):
    """ Flattens one level of a nested list.

    Parameters
    ----------
    list_to_flatten : list
        List with nested lists.

    Returns
    -------
    flattened_list : str
        Input list flattened by one level.
    """

    flattened_list = list(itertools.chain(*list_to_flatten))

    return flattened_list


def main(query_schema, subject_schema, output_path, blast_score_ratio,
         cpu_cores):

    # create output directory
    if os.path.isdir(output_path) is True:
        sys.exit('Output directory exists. Please provide '
                 'a path to a new directory that will be '
                 'created to store the files created by '
                 'the process.')
    else:
        os.mkdir(output_path)

    # import representative sequences in query schema
    rep_dir = os.path.join(query_schema, 'short')
    rep_files = [os.path.join(rep_dir, f)
                 for f in os.listdir(rep_dir)
                 if f.endswith('.fasta') is True]

    # get representative sequences from query schema
    query_ids = [os.path.basename(f).split('_')[0] for f in rep_files]
    query_reps = []
    for f in rep_files:
        locus_id = os.path.basename(f).split('_short')[0]
        records = SeqIO.parse(f, 'fasta')
        # only get the forst representative allele
        rec = next(records, None)
        seqid = rec.id
        allele_id = seqid.split('_')[-1]
        short_seqid = '{0}_{1}'.format(locus_id, allele_id)
        prot = translate_sequence(str(rec.seq), 11)
        sequence = '>{0}\n{1}'.format(short_seqid, prot)
        query_reps.append(sequence)

    # save query reps into same file
    query_prot_file = os.path.join(output_path, 'query_prots.fasta')
    with open(query_prot_file, 'w') as op:
        op.write('\n'.join(query_reps))

    # create BLAST db with query sequences and get self scores
    query_blastdb_path = os.path.join(output_path, 'query_blastdb')
    make_blast_db(query_prot_file, query_blastdb_path, 'prot')

    # determine self raw score for representative sequences
    self_blast_out = os.path.join(output_path, 'self_results.tsv')
    # max_targets has to be greater than 1
    # for some cases, the first alignment that is reported is
    # not the self-alignment
    run_blast('blastp', query_blastdb_path, query_prot_file, self_blast_out,
              max_hsps=1, threads=cpu_cores, ids_file=None, max_targets=5)

    self_blast_results = read_tabular(self_blast_out)
    self_blast_results = {r[0].split('_')[0]: r[2]
                          for r in self_blast_results
                          if r[0] == r[1]}

    # translate subject sequences
    subject_files = [os.path.join(subject_schema, f)
                     for f in os.listdir(subject_schema)
                     if f.endswith('.fasta') is True]

    subject_prots_file = os.path.join(output_path, 'subject_prots.fasta')
    ids = {}
    start = 1
    for file in subject_files:
        records = [(rec.id, str(rec.seq))
                   for rec in SeqIO.parse(file, 'fasta')]
        for rec in records:
            ids[rec[0]] = start
            start += 1
        sequences = ['>{0}\n{1}'.format(ids[rec[0]],
                                        translate_sequence(rec[1], 11))
                     for rec in records]
        with open(subject_prots_file, 'a') as sf:
            sf.write('\n'.join(sequences)+'\n')

    # create BLASTdb with subject sequences
    blastdb_path = os.path.join(output_path, 'subject_blastdb')
    make_blast_db(subject_prots_file, blastdb_path, 'prot')

    # BLASTp old seqs against new seqs
    blast_out = os.path.join(output_path, 'results.tsv')
    run_blast('blastp', blastdb_path, query_prot_file, blast_out,
              max_hsps=1, threads=cpu_cores, ids_file=None, max_targets=10)

    # import BLAST results
    blast_results = read_tabular(blast_out)

    ids_rev = {v: k for k, v in ids.items()}

    # determine BSR values
    bsr_values = {}
    multiple_matches = {}
    for m in blast_results:
        query = m[0].split('_')[0]
        subject = ids_rev[int(m[1])]
        score = m[-1]
        bsr = float(score) / float(self_blast_results[query])
        if query in bsr_values:
            if bsr > bsr_values[query][1]:
                bsr_values[query] = [subject, bsr]
            if bsr > blast_score_ratio:
                multiple_matches[query].append([subject, bsr])
        elif query not in bsr_values and bsr > blast_score_ratio:
            bsr_values[query] = [subject, bsr]
            multiple_matches[query] = [[subject, bsr]]

    # keep only queries with multiple matches
    multiple = []
    for k, v in multiple_matches.items():
        loci = [e[0].split('_')[0] for e in v]
        if len(set(loci)) > 1:
            matches = ['{0}\t{1}\t{2}'.format(k, e[0], e[1]) for e in v]
            multiple.extend(matches)

    multiple_lines = '\n'.join(multiple)

    multiple_file = os.path.join(output_path, 'multiple_matches.tsv')
    with open(multiple_file, 'w') as mh:
        mh.write(multiple_lines+'\n')

    # save matches between schemas loci
    matches = ['{0}\t{1}\t{2}'.format(k, v[0].split('_')[0], v[1])
               for k, v in bsr_values.items()]
    matches_lines = '\n'.join(matches)
    matches_file = os.path.join(output_path, 'matches.tsv')
    with open(matches_file, 'w') as mf:
        mf.write(matches_lines+'\n')

    # determine identifiers that had no match
    no_match = [i for i in self_blast_results if i not in bsr_values]
    no_match_lines = '\n'.join(no_match)
    no_match_file = os.path.join(output_path, 'no_match.txt')
    with open(no_match_file, 'w') as nm:
        nm.write(no_match_lines+'\n')


def parse_arguments():

    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument('-q', type=str, required=True,
                        dest='query_schema',
                        help='Path to the query schema directory.'
                             'This schema will be matched against '
                             'the subject schema.')

    parser.add_argument('-s', type=str, required=True,
                        dest='subject_schema',
                        help='Path to que subject schema directory.')

    parser.add_argument('-o', type=str, required=True,
                        dest='output_path',
                        help='Path to the output directory.')

    parser.add_argument('--bsr', type=float, required=False,
                        default=0.6, dest='blast_score_ratio',
                        help='Minimum BSR value to consider aligned '
                             'alleles as alleles for the same locus.')

    parser.add_argument('--cpu', type=int, required=False,
                        default=1, dest='cpu_cores',
                        help='Number of CPU cores to pass to BLAST.')

    args = parser.parse_args()

    return args


if __name__ == '__main__':

    args = parse_arguments()

    main(**vars(args))
