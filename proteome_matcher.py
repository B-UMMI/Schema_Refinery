#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
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

        Args:
            string (str): string to be reversed.

        Returns:
            revstr (str): reverse of input string.
    """

    revstr = string[::-1]

    return revstr


def reverse_complement(dna_sequence):
    """ Determines the reverse complement of given DNA strand.

        Args:
            dna_sequence (str): string representing a DNA sequence.

        Returns:
            reverse_complement_strand (str): the reverse complement
            of the DNA sequence, without lowercase letters.
    """

    base_complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A',
                       'a': 'T', 'c': 'G', 'g': 'C', 't': 'A'}

    # convert string into list with each character as a separate element
    bases = list(dna_sequence)

    # determine complement strand
    complement_bases = []
    for base in bases:
        if base in base_complement:
            complement_bases.append(base_complement[base])
        else:
            complement_bases.append(base.upper())

    complement_strand = ''.join(complement_bases)

    # reverse strand
    reverse_complement_strand = reverse_str(complement_strand)

    return reverse_complement_strand


def translate_sequence(dna_str, table_id):
    """ Translate a DNA sequence using the BioPython package.

        Args:
            dna_str (str): DNA sequence as string type.
            table_id (int): translation table identifier.

        Returns:
            protseq (str): protein sequence created by translating
            the input DNA sequence.
    """

    myseq_obj = Seq(dna_str)
    # sequences must be a complete and valid CDS
    protseq = Seq.translate(myseq_obj, table=table_id, cds=True)

    return protseq


def make_blast_db(makeblastdb_path, input_fasta, output_path, db_type,
                  ignore=None):
    """ Creates a BLAST database.
        Parameters
        ----------
        input_fasta : str
            Path to the FASTA file that contains the sequences
            that should be added to the BLAST database.
        output_path : str
            Path to the directory where the database files
            will be created. Database files will have names
            with the path's basemane.
        db_type : str
            Type of the database, nucleotide (nuc) or
            protein (prot).
        Returns
        -------
        Creates a BLAST database with the input sequences.
    """

    blastdb_cmd = [makeblastdb_path, '-in', input_fasta, '-out', output_path,
                   '-parse_seqids', '-dbtype', db_type]

    makedb_cmd = subprocess.Popen(blastdb_cmd,
                                  stdout=subprocess.PIPE,
                                  stderr=subprocess.PIPE)

    stderr = makedb_cmd.stderr.readlines()

    return stderr


def run_blast(blast_path, blast_db, fasta_file, blast_output,
              max_hsps=1, threads=1, ids_file=None, blast_task=None,
              max_targets=None, ignore=None):
    """ Execute BLAST to align sequences in a FASTA file
        against a BLAST database.
        Parameters
        ----------
        blast_path : str
            Path to BLAST executables.
        blast_db : str
            Path to the BLAST database.
        fasta_file : str
            Path to the FASTA file with sequences to
            align against the BLAST database.
        blast_output : str
            Path to the file that will be created to
            store BLAST results.
        max_hsps : int
            Maximum number of High Scoring Pairs per
            pair of aligned sequences.
        threads : int
            Number of threads/cores used to run BLAST.
        ids_file : str
            Path to a file with sequence identifiers,
            one per line. Sequences will only be aligned
            to the sequences in the BLAST database that
            have any of the identifiers in this file.
        blast_task : str
            Type of BLAST task.
        max_targets : int
            Maximum number of target of subject sequences
            to align against.
        Returns
        -------
        stderr : str
            String with errors raised during BLAST execution.
    """

    blast_args = [blast_path, '-db', blast_db, '-query', fasta_file,
                  '-out', blast_output, '-outfmt', '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore score',
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


#schema_dir = '/home/rfm/Desktop/rfm/Lab_Analyses/GAS_PrepExternalSchema/Annotations/schema_seed_gas_dys'
#references_dir = '/home/rfm/Desktop/rfm/Lab_Analyses/GAS_PrepExternalSchema/Annotations/splitted_proteomes'
#output_dir = '/home/rfm/Desktop/rfm/Lab_Analyses/GAS_PrepExternalSchema/Annotations/proteome_matching'
#blast_threads = 6
def main(schema_dir, references_dir, output_dir, blast_threads):

    if os.path.isdir(output_dir) is False:
        os.mkdir(output_dir)

    reps_dir = os.path.join(schema_dir, 'short')
    rep_files = [os.path.join(reps_dir, f) for f in os.listdir(reps_dir) if f.endswith('.fasta') is True]

    # get all representative sequences into same file
    reps = []
    for f in rep_files:
        for rec in SeqIO.parse(f, 'fasta'):
            prot = translate_sequence(str(rec.seq), 11)
            sequence = '>{0}\n{1}'.format(rec.id, prot)
            reps.append(sequence)

    # save new reps into same file
    prot_file = os.path.join(output_dir, 'reps_prots.fasta')
    with open(prot_file, 'w') as pinfile:
        pinfile.write('\n'.join(reps))

    # BLASTp TrEMBL and Swiss-Prot records

    # create TrEMBL BLASTdb
    tr_file = os.path.join(references_dir, 'trembl_prots.fasta')
    tr_blastdb_path = os.path.join(output_dir, 'tr_db')
    tr_blastdb_stderr = make_blast_db('makeblastdb', tr_file, tr_blastdb_path, 'prot')
    tr_blastout = os.path.join(output_dir, 'tr_blastout.tsv')
    tr_blast_stderr = run_blast('blastp', tr_blastdb_path, prot_file,
                                tr_blastout, max_hsps=1, threads=blast_threads,
                                ids_file=None, blast_task=None, max_targets=1)

    # create Swiss-Prot BASLTdb
    sp_file = os.path.join(references_dir, 'sp_prots.fasta')
    sp_blastdb_path = os.path.join(output_dir, 'sp_db')
    sp_blastdb_stderr = make_blast_db('makeblastdb', sp_file, sp_blastdb_path, 'prot')
    sp_blastout = os.path.join(output_dir, 'sp_blastout.tsv')
    sp_blast_stderr = run_blast('blastp', sp_blastdb_path, prot_file,
                                sp_blastout, max_hsps=1, threads=blast_threads,
                                ids_file=None, blast_task=None, max_targets=1)

    # self BLASTp to get reps self Raw Scores
    reps_blastdb_path = os.path.join(output_dir, os.path.join(output_dir, 'reps_db'))
    reps_blastdb_stderr = make_blast_db('makeblastdb', prot_file, reps_blastdb_path, 'prot')
    reps_self_blastout = os.path.join(output_dir, 'reps_self_blastout.tsv')
    reps_blast_stderr = run_blast('blastp', reps_blastdb_path, prot_file,
                                  reps_self_blastout, max_hsps=1, threads=blast_threads,
                                  ids_file=None, blast_task=None, max_targets=1)

    # import self_results
    with open(reps_self_blastout, 'r') as at:
        reps_self_results = list(csv.reader(at, delimiter='\t'))

    reps_self_scores = {l[0]: l[-1] for l in reps_self_results}

    # import Swiss-Prot and TrEMBL records descriptions
    with open(os.path.join(references_dir, 'descriptions'), 'rb') as dinfile:
        descriptions = pickle.load(dinfile)

    # get TrEMBL and Swiss-Prot results and choose only highest representative hit for each gene    
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
    tr_header = 'TrEMBL_ID\tTrEMBL_BSR\tTrEMBL_LNAME\tTrEMBL_SNAME'
    tr_annotations = os.path.join(output_dir, 'tr_annotations.tsv')
    with open(tr_annotations, 'w') as trout:
        tr_outlines = [tr_header] + ['{0}\t{1}\t{2}\t{3}'.format(*v) for k, v in tr_selected.items()]
        tr_outtext = '\n'.join(tr_outlines)
        trout.write(tr_outtext)

    sp_header = 'SwissProt_ID\tSwissProt_BSR\tSwissProt_LNAME\tSwissProt_SNAME'
    sp_annotations = os.path.join(output_dir, 'sp_annotations.tsv')
    with open(sp_annotations, 'w') as spout:
        sp_outlines = [sp_header] + ['{0}\t{1}\t{2}\t{3}'.format(*v) for k, v in sp_selected.items()]
        sp_outtext = '\n'.join(sp_outlines)
        spout.write(sp_outtext)


def parse_arguments():

    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument('-s', type=str, required=True,
                        dest='schema_dir',
                        help='Path to the schema\'s directory.')

    parser.add_argument('-r', type=str, required=True,
                        dest='references_dir',
                        help='Path to the directory with Swiss-Prot and'
                             'and file with descriptions.')

    parser.add_argument('-o', type=str, required=True,
                        dest='output_dir',
                        help='Path to output directory.')

    parser.add_argument('--bt', type=int, required=False,
                        dest='blast_threads',
                        help='Number of threads to pass to BLAST.')

    args = parser.parse_args()

    return args


if __name__ == '__main__':

    args = parse_arguments()

    main(**vars(args))
