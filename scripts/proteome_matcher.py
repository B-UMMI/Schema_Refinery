#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Purpose
-------
This script translates loci representative alleles and
aligns them against Swiss-Prot and TrEMBL records to
select annotation terms based on the BSR computed for
each alignment.

Code documentation
------------------
"""


import os
import csv
import pickle
import argparse
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


def main(schema_directory, records_directory, output_directory, cpu_cores):

    if os.path.isdir(output_directory) is False:
        os.mkdir(output_directory)

    reps_dir = os.path.join(schema_directory, 'short')
    rep_files = [os.path.join(reps_dir, f)
                 for f in os.listdir(reps_dir)
                 if f.endswith('.fasta') is True]

    # get all representative sequences into same file
    reps = []
    for f in rep_files:
        locus_id = os.path.basename(f).split('_short')[0]
        for rec in SeqIO.parse(f, 'fasta'):
            seqid = rec.id
            allele_id = seqid.split('_')[-1]
            short_seqid = '{0}_{1}'.format(locus_id, allele_id)
            prot = translate_sequence(str(rec.seq), 11)
            sequence = '>{0}\n{1}'.format(short_seqid, prot)
            reps.append(sequence)

    # save new reps into same file
    prot_file = os.path.join(output_directory, 'reps_prots.fasta')
    with open(prot_file, 'w') as pinfile:
        pinfile.write('\n'.join(reps)+'\n')

    # BLASTp TrEMBL and Swiss-Prot records

    # create TrEMBL BLASTdb
    tr_file = os.path.join(records_directory, 'trembl_prots.fasta')
    tr_blastdb_path = os.path.join(output_directory, 'tr_db')
    tr_blastdb_stderr = make_blast_db(tr_file, tr_blastdb_path, 'prot')
    tr_blastout = os.path.join(output_directory, 'tr_blastout.tsv')
    tr_blast_stderr = run_blast('blastp', tr_blastdb_path, prot_file,
                                tr_blastout, max_hsps=1, threads=cpu_cores,
                                ids_file=None, blast_task=None, max_targets=1)

    # create Swiss-Prot BASLTdb
    sp_file = os.path.join(records_directory, 'sp_prots.fasta')
    sp_blastdb_path = os.path.join(output_directory, 'sp_db')
    sp_blastdb_stderr = make_blast_db(sp_file, sp_blastdb_path, 'prot')
    sp_blastout = os.path.join(output_directory, 'sp_blastout.tsv')
    sp_blast_stderr = run_blast('blastp', sp_blastdb_path, prot_file,
                                sp_blastout, max_hsps=1, threads=cpu_cores,
                                ids_file=None, blast_task=None, max_targets=1)

    # self BLASTp to get reps self Raw Scores
    reps_blastdb_path = os.path.join(output_directory, os.path.join(output_directory, 'reps_db'))
    reps_blastdb_stderr = make_blast_db(prot_file, reps_blastdb_path, 'prot')
    reps_self_blastout = os.path.join(output_directory, 'reps_self_blastout.tsv')
    reps_blast_stderr = run_blast('blastp', reps_blastdb_path, prot_file,
                                  reps_self_blastout, max_hsps=1, threads=cpu_cores,
                                  ids_file=None, blast_task=None, max_targets=10)

    # import self_results
    with open(reps_self_blastout, 'r') as at:
        reps_self_results = list(csv.reader(at, delimiter='\t'))

    reps_self_scores = {l[0]: l[-1]
                        for l in reps_self_results
                        if l[0] == l[1]}

    # import Swiss-Prot and TrEMBL records descriptions
    with open(os.path.join(records_directory, 'descriptions'), 'rb') as dinfile:
        descriptions = pickle.load(dinfile)

    # get TrEMBL and Swiss-Prot results and choose only
    # highest representative hit for each gene
    with open(tr_blastout, 'r') as at:
        tr_lines = list(csv.reader(at, delimiter='\t'))

    # TrEMBL
    tr_results = {}
    for l in tr_lines:
        query = l[0]
        subject = l[1]
        score = l[-1]
        self_score = reps_self_scores[query]
        bsr = float(score) / float(self_score)
        locus = query.split('_')[0]
        tr_results.setdefault(locus, []).append([subject, bsr])

    # select only best hit
    for k, v in tr_results.items():
        if len(v) > 1:
            sorted_v = sorted(v, key= lambda x: float(x[1]), reverse=True)
            tr_results[k] = sorted_v[0]
        else:
            tr_results[k] = v[0]

    # get short and long names from description
    tr_selected = {k: v for k, v in tr_results.items() if v[1] >= 0.6}
    for k, v in tr_selected.items():
        desc = descriptions[v[0]]
        lname = desc.split(v[0]+' ')[1].split(' OS=')[0]
        sname = desc.split('GN=')[1].split(' PE=')[0]
        tr_selected[k] = v + [lname, sname]

    # Swiss-Prot
    with open(sp_blastout, 'r') as at:
        sp_lines = list(csv.reader(at, delimiter='\t'))

    # Swiss-Prot
    sp_results = {}
    for l in sp_lines:
        query = l[0]
        subject = l[1]
        score = l[-1]
        self_score = reps_self_scores[query]
        bsr = float(score) / float(self_score)
        locus = query.split('_')[0]
        sp_results.setdefault(locus, []).append([subject, bsr])

    for k, v in sp_results.items():
        if len(v) > 1:
            sorted_v = sorted(v, key=lambda x: float(x[1]), reverse=True)
            sp_results[k] = sorted_v[0]
        else:
            sp_results[k] = v[0]

    sp_selected = {k: v for k, v in sp_results.items() if v[1] >= 0.6}
    for k, v in sp_selected.items():
        desc = descriptions[v[0]]
        lname = desc.split(v[0]+' ')[1].split(' OS=')[0]
        sname = desc.split('GN=')[1].split(' PE=')[0]
        sp_selected[k] = v + [lname, sname]

    # save results
    tr_header = 'Locus_ID\tTrEMBL_ID\tTrEMBL_BSR\tTrEMBL_LNAME\tTrEMBL_SNAME'
    tr_annotations = os.path.join(output_directory, 'tr_annotations.tsv')
    with open(tr_annotations, 'w') as trout:
        tr_outlines = [tr_header] + ['{0}\t{1}\t{2}\t{3}\t{4}'.format(k,*v) for k, v in tr_selected.items()]
        tr_outtext = '\n'.join(tr_outlines)
        trout.write(tr_outtext+'\n')

    print('TrEMBL annotations available at {0}'.format(tr_annotations))

    sp_header = 'Locus_ID\tSwissProt_ID\tSwissProt_BSR\tSwissProt_LNAME\tSwissProt_SNAME'
    sp_annotations = os.path.join(output_directory, 'sp_annotations.tsv')
    with open(sp_annotations, 'w') as spout:
        sp_outlines = [sp_header] + ['{0}\t{1}\t{2}\t{3}\t{4}'.format(k,*v) for k, v in sp_selected.items()]
        sp_outtext = '\n'.join(sp_outlines)
        spout.write(sp_outtext+'\n')

    print('SwissProt annotations available at {0}'.format(sp_annotations))


def parse_arguments():

    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument('-s', '--schema-directory', type=str,
                        required=True, dest='schema_directory',
                        help='Path to the schema\'s directory.')

    parser.add_argument('-p', '--records-directory', type=str, required=True,
                        dest='records_directory',
                        help='Path to the directory with the Fasta files '
                             'with Swiss-Prot and TrEMBL records. It must '
                             'also contain a file with record descriptions.')

    parser.add_argument('-o', '--output-directory', type=str,
                        required=True, dest='output_directory',
                        help='Path to the output directory.')

    parser.add_argument('-cpu', '--cpu-cores', type=int, required=False,
                        dest='cpu_cores',
                        default=1,
                        help='Number of CPU cores to pass to BLAST.')

    args = parser.parse_args()

    return args


if __name__ == '__main__':

    args = parse_arguments()
    main(**vars(args))
