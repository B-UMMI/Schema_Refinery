#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

"""


import os
import csv
import argparse
import subprocess

import networkx as nx
from Bio import SeqIO
from Bio.Seq import Seq
from pyvis.network import Network


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


#networkx_graph = G
#reps = representatives
def graph_html(networkx_graph, reps, output_html, loci_ids, title):
    """
    """

    nt = Network('720px', '1280px')

    # attribute colors to each locus. Hardcoded!
    colors = [('#9ecae1', '#2171b5'), ('#a1d99b', '#238b45')]
    colors_map = {locus: colors[i] for i, locus in enumerate(loci_ids)}

    # populates the nodes and edges data structures
    nt.from_nx(networkx_graph)

    for node in nt.nodes:
        node_id = node['id']
        node_locus, allele_id = node_id.split('_')
        node['label'] = allele_id
        if node_id not in reps:
            node['color'] = colors_map[node_locus][0]
        else:
            node['color'] = colors_map[node_locus][1]
            node['size'] = 12

#    for edge in nt.edges:
#        rounded_weight = round(edge['weight'], 2)
#        edge['title'] = str(rounded_weight)

    nt.show_buttons(filter_=['physics', 'manipulation'])

    # change default title
    nt.heading = title
    # enable multiselection of nodes and edges
    nt.options.interaction.multiselect = True
    # enable edge highlighting during hover
    nt.options.interaction.hover = True
    # change deafults for BarnesHut
#    nt.barnes_hut(gravity=-2000, central_gravity=0.3, spring_length=95,
#                  spring_strength=0.04, damping=0.25, overlap=0.5)

    nt.show(output_html)


# add way to apply cutoff and output groups of alleles
loci_ids = ['GCF-002236855-protein1020']
schema_dir = '/home/rfm/Desktop/rfm/Lab_Analyses/GAS_PrepExternalSchema/wgMLST_schema/spyogenes_schema/solve_problematic_2/spyogenes_schema_processed'
output_dir = '/home/rfm/Desktop/rfm/Lab_Analyses/GAS_PrepExternalSchema/wgMLST_schema/spyogenes_schema/solve_problematic_2/scl_analysis/test_graph'
blast_score_ratio = 0.6
title = 'sclB'
def main(loci_ids, schema_dir, output_dir, blast_score_ratio,
         title, frequency_table):

    if os.path.isdir(output_dir) is False:
        os.mkdir(output_dir)
    
    loci_files = {locus: os.path.join(schema_dir, locus+'.fasta') for locus in loci_ids}

    # create file with short ids
    records = {k: [] for k in loci_files}
    for locus, file in loci_files.items():
        current_records = [('{0}_{1}'.format(locus, (rec.id).split('_')[-1]), str((rec.translate(table=11, cds=True)).seq))
                           for rec in SeqIO.parse(file, 'fasta')]
        records[locus] = current_records

    # get locus rep IDs
    loci_rep_files = {locus: os.path.join(schema_dir, 'short', locus+'_short.fasta') for locus in loci_ids}
    representatives = []
    for locus, file in loci_rep_files.items():
        current_representatives = ['{0}_{1}'.format(locus, (rec.id).split('_')[-1]) for rec in SeqIO.parse(file, 'fasta')]
        representatives.extend(current_representatives)

    # save records with short ids
    seqs_file = os.path.join(output_dir, 'locus_seqs.fasta')
    for locus, recs in records.items():
        records_str = '\n'.join(['>{0}\n{1}'.format(*rec) for rec in recs]) + '\n'
        with open(seqs_file, 'a') as outfile:
            outfile.write(records_str)

    # create BLASTdb with locus proteins
    output_path = os.path.join(output_dir, 'locus_db')
    make_blast_db('makeblastdb', seqs_file, output_path, 'prot')

    # BLAST all against all
    blast_output = os.path.join(output_dir, 'blastout.tsv')
    run_blast('blastp', output_path, seqs_file, blast_output,
              max_hsps=1, threads=6, ids_file=None, blast_task=None,
              max_targets=10, ignore=None)

    # read BLASTp results
    with open(blast_output, 'r') as infile:
        blast_results = list(csv.reader(infile, delimiter='\t'))

    # self scores
    self_scores = {l[0]: l[-1] for l in blast_results if l[0] == l[1]}

    # exclude self hits from blast results
    blast_results_filtered = [l for l in blast_results if l[0] != l[1]]

    # exclude based on evalue
    blast_results_filtered = [l for l in blast_results_filtered if float(l[10]) <= 10**-13]

    # compute and add BSR
    for i in range(len(blast_results_filtered)):
        query = blast_results_filtered[i][0]
        score = blast_results_filtered[i][-1]
        bsr = float(score) / float(self_scores[query])
        blast_results_filtered[i].append(bsr)

    # exclude based on BSR
    blast_results_filtered = [l for l in blast_results_filtered if float(l[-1]) >= blast_score_ratio]

    # log transform frequency values
    if frequency_table is not None:
        with open(frequency_table, 'r') as infile:
            frequencies = list(csv.reader(infile, delimiter='\t'))

    # instantiate graph
    G = nx.Graph()

    # add edges with weight
    edges = [(l[0], l[1], l[-1]) for l in blast_results_filtered]
    G.add_weighted_edges_from(edges)

    # create HTML with grpah viz
    output_html = os.path.join(output_dir, '{0}_graph.html'.format(title))
    graph_html(G, representatives, output_html, loci_ids, title)


# add option to scan for motifs in alleles and color nodes based on that information
def parse_arguments():

    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument('-l', type=str, required=True,
                        nargs='+', dest='loci_ids',
                        help='')

    parser.add_argument('-s', type=str, required=True,
                        dest='schema_dir',
                        help='')

    parser.add_argument('-o', type=str, required=True,
                        dest='output_dir',
                        help='')

    parser.add_argument('--bsr', type=float, required=False,
                        default=0.6,
                        dest='blast_score_ratio',
                        help='')
    
    parser.add_argument('--t', type=str, required=False,
                        default='myGraph',
                        dest='title',
                        help='')

    parser.add_argument('--f', type=str, required=False,
                        dest='frequency_table',
                        help='')

    args = parser.parse_args()

    return args


if __name__ == '__main__':

    args = parse_arguments()

    main(**vars(args))
