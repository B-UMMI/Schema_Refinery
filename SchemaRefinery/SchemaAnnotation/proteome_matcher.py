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
    from utils.blast_functions import make_blast_db, run_blast
    from utils.sequence_functions import translate_sequence

except ModuleNotFoundError:
    from SchemaRefinery.utils.blast_functions import make_blast_db, run_blast
    from SchemaRefinery.utils.sequence_functions import translate_sequence


# proteome_splitter_files_list has the TrEMBL file in index 0,
# Swiss-Prot in index 1 and descriptions file in index 2
def proteome_matcher(schema_directory:str, proteome_splitter_files_list:list, output_directory:str, cpu_cores:int):

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
    tr_file = proteome_splitter_files_list[0]
    tr_blastdb_path = os.path.join(output_directory, 'tr_db')
    tr_blastdb_stderr = make_blast_db(tr_file, tr_blastdb_path, 'prot')
    tr_blastout = os.path.join(output_directory, 'tr_blastout.tsv')
    tr_blast_stderr = run_blast('blastp', tr_blastdb_path, prot_file,
                                tr_blastout, max_hsps=1, threads=cpu_cores,
                                ids_file=None, blast_task=None, max_targets=1)

    # create Swiss-Prot BASLTdb
    sp_file = proteome_splitter_files_list[1]
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
    with open(proteome_splitter_files_list[2], 'rb') as dinfile:
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

    return tr_annotations, sp_annotations
