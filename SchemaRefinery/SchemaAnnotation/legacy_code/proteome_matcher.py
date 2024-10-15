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

from Bio import SeqIO

try:
    from utils import blast_functions as bf
    from utils.sequence_functions import translate_sequence
except ModuleNotFoundError:
    from SchemaRefinery.utils import blast_functions as bf
    from SchemaRefinery.utils.sequence_functions import translate_sequence


# proteome_files: TrEMBL FASTA, Swiss-Prot FASTA and descriptions pickle
def proteome_matcher(schema_directory: str, proteome_files: list,
                     output_directory: str, cpu_cores: int,
                     bsr: float):

    reps_dir = os.path.join(schema_directory, 'short')
    rep_files = [os.path.join(reps_dir, f)
                 for f in os.listdir(reps_dir)
                 if f.endswith('.fasta')]

    # Translate representative alleles and save to single FASTA file
    translated_reps = []
    for f in rep_files:
        locus_id = os.path.basename(f).split('_short')[0]
        for rec in SeqIO.parse(f, 'fasta'):
            seqid = rec.id
            allele_id = seqid.split('_')[-1]
            short_seqid = '{0}_{1}'.format(locus_id, allele_id)
            prot = translate_sequence(str(rec.seq), 11)
            prot_record = '>{0}\n{1}'.format(short_seqid, prot)
            translated_reps.append(prot_record)

    ouput_text = '\n'.join(translated_reps)
    prot_file = os.path.join(output_directory, 'reps_prots.fasta')
    with open(prot_file, 'w') as pinfile:
        pinfile.write(ouput_text+'\n')

    # BLASTp TrEMBL and Swiss-Prot records
    # Create TrEMBL BLASTdb
    tr_file = proteome_files[0]
    tr_blastdb_path = os.path.join(output_directory, 'tr_BLASTdb')
    tr_blastdb_stderr = bf.make_blast_db(tr_file, tr_blastdb_path, 'prot')
    tr_blastout = os.path.join(output_directory, 'tr_blastout.tsv')
    tr_blast_stderr = bf.run_blast('blastp', tr_blastdb_path, prot_file,
                                   tr_blastout, max_hsps=1, threads=cpu_cores,
                                   max_targets=1)

    # Create Swiss-Prot BLASTdb
    sp_file = proteome_files[1]
    sp_blastdb_path = os.path.join(output_directory, 'sp_db')
    sp_blastdb_stderr = bf.make_blast_db(sp_file, sp_blastdb_path, 'prot')
    sp_blastout = os.path.join(output_directory, 'sp_blastout.tsv')
    sp_blast_stderr = bf.run_blast('blastp', sp_blastdb_path, prot_file,
                                   sp_blastout, max_hsps=1, threads=cpu_cores,
                                   max_targets=1)

    # Self BLASTp to get representatives self scores
    reps_blastdb_path = os.path.join(output_directory, 'reps_db')
    reps_blastdb_stderr = bf.make_blast_db(prot_file, reps_blastdb_path, 'prot')
    reps_self_blastout = os.path.join(output_directory, 'reps_self_blastout.tsv')
    reps_blast_stderr = bf.run_blast('blastp', reps_blastdb_path, prot_file,
                                     reps_self_blastout, max_hsps=1, threads=cpu_cores,
                                     max_targets=10)

    # Import self_results
    with open(reps_self_blastout, 'r') as at:
        results = list(csv.reader(at, delimiter='\t'))

    self_scores = {line[0]: line[-1] for line in results if line[0] == line[1]}

    # Import Swiss-Prot and TrEMBL records descriptions
    with open(proteome_files[2], 'rb') as dinfile:
        descriptions = pickle.load(dinfile)

    # Get TrEMBL and Swiss-Prot results and choose only highest-score matches
    with open(tr_blastout, 'r') as at:
        tr_lines = list(csv.reader(at, delimiter='\t'))

    # TrEMBL
    tr_results = {}
    for line in tr_lines:
        query = line[0]
        subject = line[1]
        score = line[-1]
        self_score = self_scores[query]
        match_bsr = float(score) / float(self_score)
        locus = query.split('_')[0]
        tr_results.setdefault(locus, []).append([subject, match_bsr])

    # Select only best hit
    for k, v in tr_results.items():
        if len(v) > 1:
            sorted_v = sorted(v, key=lambda x: float(x[1]), reverse=True)
            tr_results[k] = sorted_v[0]
        else:
            tr_results[k] = v[0]

    # Get short and long names from description
    tr_selected = {k: v for k, v in tr_results.items() if v[1] >= bsr}
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
    for line in sp_lines:
        query = line[0]
        subject = line[1]
        score = line[-1]
        self_score = self_scores[query]
        match_bsr = float(score) / float(self_score)
        locus = query.split('_')[0]
        sp_results.setdefault(locus, []).append([subject, match_bsr])

    for k, v in sp_results.items():
        if len(v) > 1:
            sorted_v = sorted(v, key=lambda x: float(x[1]), reverse=True)
            sp_results[k] = sorted_v[0]
        else:
            sp_results[k] = v[0]

    sp_selected = {k: v for k, v in sp_results.items() if v[1] >= bsr}
    for k, v in sp_selected.items():
        desc = descriptions[v[0]]
        lname = desc.split(v[0]+' ')[1].split(' OS=')[0]
        sname = desc.split('GN=')[1].split(' PE=')[0]
        sp_selected[k] = v + [lname, sname]

    # Save results
    tr_header = 'Locus\tTrEMBL_ID\tTrEMBL_BSR\tTrEMBL_LNAME\tTrEMBL_SNAME'
    tr_annotations = os.path.join(output_directory, 'tr_annotations.tsv')
    with open(tr_annotations, 'w') as trout:
        tr_outlines = [tr_header] + ['{0}\t{1}\t{2}\t{3}\t{4}'.format(k, *v) for k, v in tr_selected.items()]
        tr_outtext = '\n'.join(tr_outlines)
        trout.write(tr_outtext+'\n')

    print('TrEMBL annotations available at {0}'.format(tr_annotations))

    sp_header = 'Locus\tSwissProt_ID\tSwissProt_BSR\tSwissProt_LNAME\tSwissProt_SNAME'
    sp_annotations = os.path.join(output_directory, 'sp_annotations.tsv')
    with open(sp_annotations, 'w') as spout:
        sp_outlines = [sp_header] + ['{0}\t{1}\t{2}\t{3}\t{4}'.format(k, *v) for k, v in sp_selected.items()]
        sp_outtext = '\n'.join(sp_outlines)
        spout.write(sp_outtext+'\n')

    print('SwissProt annotations available at {0}'.format(sp_annotations))

    return tr_annotations, sp_annotations
