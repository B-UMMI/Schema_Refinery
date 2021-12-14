#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Purpose
-------
This script creates a HTML file with a network. Nodes are
linked based on the BSR value computed for each alignment
between all translated sequences in the input Fasta file.
The input Fasta file must be a Fasta file from a chewBBACA
schema with the translated alleles.

Code documentation
------------------
"""


import os
import csv
import argparse
import subprocess

import networkx as nx
from Bio import SeqIO
from pyvis.network import Network


def make_blast_db(makeblastdb_path, input_fasta, output_path, db_type):
    """ Creates a BLAST database.

    Parameters
    ----------
    makeblastdb_path :  str
        Path to the makeblastdb executable.
    input_fasta : str
        Path to the FASTA file that contains the sequences
        that should be added to the BLAST database.
    output_path : str
        Path to the directory where the database files
        will be created. Database files will have names
        with the path's basemane.
    db_type : str
        Database type , nucleotide (nuc) or protein (prot).

    Returns
    -------
    Creates a BLAST database with the input sequences.
    """

    blastdb_cmd = [makeblastdb_path, '-in', input_fasta,
                   '-out', output_path, '-parse_seqids',
                   '-dbtype', db_type]

    makedb_cmd = subprocess.Popen(blastdb_cmd,
                                  stdout=subprocess.PIPE,
                                  stderr=subprocess.PIPE)

    stderr = makedb_cmd.stderr.readlines()

    return stderr


def run_blast(blast_path, blast_db, fasta_file, blast_output,
              max_hsps=1, threads=1, ids_file=None, blast_task=None,
              max_targets=None):
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
        Maximum number of target subject sequences
        to align against.

    Returns
    -------
    stderr : str
        String with errors raised during BLAST execution.
    """

    outfmt = ('6 qseqid sseqid pident length mismatch gapopen '
              'qstart qend sstart send evalue bitscore score')
    blast_args = [blast_path, '-db', blast_db, '-query', fasta_file,
                  '-out', blast_output, '-outfmt', outfmt,
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


def linear_conversion(value, current_min, current_max, new_min, new_max):
    """ Converts one value to a value in another linear scale.

    Parameters
    ----------
    value : int or float
        Value to convert to new scale.
    current_min : int or float
        Minimum value of the current scale.
    current_max : int or float
        Maximum value of the current scale.
    new_min : int or float
        Minimum value of the new scale.
    new_max : int or float
        Maximum value of the new scale.

    Returns
    -------
    converted_value : int or float
        Value converted to the new scale.
    """

    converted_value = ((value-current_min) / (current_max-current_min)) \
        * (new_max - new_min) + new_min

    return converted_value


def convert_frequencies(frequencies, new_min, new_max):
    """ Converts a set of values to new values in another linear scale.

    Parameters
    ----------
    frequencies : list or set
        A list or a set with all distinct frequency values.
    new_min : int
        Minimum value of the new scale.
    new_max : int
        Maximum value of the new scale.

    Returns
    -------
    converted_values : dict
        Dictionary with frequency values as keys and the value
        in the new scale as values.
    """

    max_freq = max(frequencies)
    min_freq = min(frequencies)
    converted_values = {f: linear_conversion(f, min_freq, max_freq, new_min, new_max)
                        for f in frequencies}

    return converted_values


def graph_html(networkx_graph, reps, output_html, title, frequencies):
    """ Create a HTML file with a network representation.

    Parameters
    ----------
    networkx_graph : networkx.classes.graph.Graph
        Networkx graph object with data to add to
        Pyvis object.
    reps : list
        List with the identifiers of the representative
        sequences.
    output_html : 
        Path to the output HTML file.
    title : str
        Title to display in the HTML page.
    frequencies : dict
        
    """

    # linear transform frequencies to control node size
    if frequencies is not None:
        # min and max values for node size
        min_size = 10
        max_size = 50
        freqs_values = set([int(v[0])
                            for k, v in frequencies.items()
                            if int(v[0]) != 0])

        # compute values for node size based on frequencies
        converted_values = convert_frequencies(freqs_values, min_size, max_size)

    # network/plot size
    nt = Network('720px', '1280px')

    # populates the nodes and edges data structures
    nt.from_nx(networkx_graph)

    for node in nt.nodes:
        node_id = node['id']

        if frequencies is not None:
            allele_ids = frequencies[node_id][1]
            allele_ids = allele_ids.split(',')
            node['title'] = ('freq: {0}<br>alleles: '
                             '{1}').format(frequencies[node_id][0],
                                           frequencies[node_id][1])

            if any([i in reps for i in allele_ids]) is True:
                node['color'] = '#a50f15'
            else:
                node['color'] = '#2171b5'

            prot_freq = int(frequencies[node_id][0])
            if prot_freq == 0:
                node['size'] = 1
            else:
                node['size'] = converted_values[prot_freq]
        # without frequency values
        else:
            node['color'] = '#2171b5'
            node['size'] = 10

    nt.show_buttons(filter_=['physics', 'manipulation'])

    # change default title
    nt.heading = title
    # enable multiselection of nodes and edges
    nt.options.interaction.multiselect = True
    # enable edge highlighting during hover
    nt.options.interaction.hover = True

    nt.show(output_html)


def main(input_file, output_dir, schema_dir, locus_id, blast_score_ratio,
         title, frequency_table, maximum_links):

    if os.path.isdir(output_dir) is False:
        os.mkdir(output_dir)

    # create Fasta file with short headers to avoid BLAST error
    # related with long sequence identifiers
    records = []
    for rec in SeqIO.parse(input_file, 'fasta'):
        short_id = '{0}_{1}'.format(locus_id, (rec.id).split('_')[-1])
        records.append('>{0}\n{1}'.format(short_id, str(rec.seq)))

    input_file_basename = os.path.basename(input_file)
    input_file_basename = input_file_basename.split('.fasta')[0] + '_shortIDs.fasta'
    input_file = os.path.join(output_dir, input_file_basename)
    with open(input_file, 'w') as outfile:
        text = '\n'.join(records)
        outfile.write(text+'\n')

    # create BLASTdb with locus proteins
    output_path = os.path.join(output_dir, 'locus_db')
    stderr = make_blast_db('makeblastdb', input_file, output_path, 'prot')
    if len(stderr) > 0:
        print(stderr)

    # BLAST all against all
    blast_output = os.path.join(output_dir, 'blastout.tsv')
    # do not set max_targets
    # self-alignment might not be reported
    stderr = run_blast('blastp', output_path, input_file,
                       blast_output, max_hsps=1, threads=6,
                       ids_file=None, blast_task=None)
    if len(stderr) > 0:
        print(stderr)

    # read BLASTp results
    with open(blast_output, 'r') as infile:
        blast_results = list(csv.reader(infile, delimiter='\t'))

    # self scores
    self_scores = {l[0]: l[-1]
                   for l in blast_results
                   if l[0] == l[1]}

    # exclude self hits from blast results
    blast_results_filtered = [l
                              for l in blast_results
                              if l[0] != l[1]]

    # exclude based on evalue
    blast_results_filtered = [l
                              for l in blast_results_filtered
                              if float(l[10]) <= 10**-13]

    # compute and add BSR
    for i in range(len(blast_results_filtered)):
        query = blast_results_filtered[i][0]
        score = blast_results_filtered[i][-1]
        bsr = float(score) / float(self_scores[query])
        blast_results_filtered[i].append(bsr)

    # exclude based on BSR
    blast_results_filtered = [l
                              for l in blast_results_filtered
                              if float(l[-1]) >= blast_score_ratio]

    # select highest scoring links for each protein
    selected_results = {}
    for r in blast_results_filtered:
        if r[0] in selected_results:
            if len(selected_results[r[0]]) < maximum_links:
                selected_results[r[0]].append(r)
        else:
            selected_results[r[0]] = [r]
            
    blast_results_filtered = []
    for k, v in selected_results.items():
        blast_results_filtered.extend(v)

    # get IDs of representatives
    rep_file = os.path.join(schema_dir, 'short', locus_id+'_short.fasta')
    representative_ids = [(rec.id).split('_')[-1]
                          for rec in SeqIO.parse(rep_file, 'fasta')]

    # read and log transform frequency values
    # not working!
    # frequencies_dict = None
    # if frequency_table is not None:
    #     with open(frequency_table, 'r') as infile:
    #         frequencies = list(csv.reader(infile, delimiter='\t'))[1:]
    #         frequencies_dict = {l[0]: l[1:] for l in frequencies}

    # instantiate graph
    G = nx.Graph()

    # add edges with weight
    edges = [(l[0], l[1], l[-1]) for l in blast_results_filtered]
    G.add_weighted_edges_from(edges)

    # create HTML with grpah viz
    output_html = os.path.join(output_dir, '{0}_graph.html'.format(title))
    graph_html(G, representative_ids, output_html, title, frequencies_dict)


def parse_arguments():

    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument('-i', type=str, required=True,
                        dest='input_file',
                        help='Fasta file with locus translated sequences.')

    parser.add_argument('-o', type=str, required=True,
                        dest='output_dir',
                        help='Path to output directory where output files '
                             'will be saved.')

    parser.add_argument('-s', type=str, required=True,
                        dest='schema_dir',
                        help='Path to the schema\'s directory.' )

    parser.add_argument('-l', type=str, required=True,
                        dest='locus_id',
                        help='The ID of the locus.')

    parser.add_argument('--bsr', type=float, required=False,
                        default=0.6,
                        dest='blast_score_ratio',
                        help='Minimum BLAST Score Ratio value. '
                             'Smaller values will not be shown as '
                             'edges in the graph.')

    parser.add_argument('--t', type=str, required=False,
                        default='myGraph',
                        dest='title',
                        help='Title to display in the HTML page.')

    # currently not working
    parser.add_argument('--f', type=str, required=False,
                        dest='frequency_table',
                        help='TSV file with protein and allele '
                             'frequencies. The first column has '
                             'the ID of the proteins in the input '
                             'Fasta file. The second column has '
                             'the frequency of each protein based '
                             'on the frequency of alleles that code '
                             'for the protein. The last column has '
                             'the identifiers of the alleles that '
                             'code for the protein (separated by commas).')

    parser.add_argument('-m', '--maximum-links', type=int,
                        required=False, default=5,
                        dest='maximum_links',
                        help='Maximum number high-scoring alignments'
                             ' selected per protein.')

    args = parser.parse_args()

    return args


if __name__ == '__main__':

    args = parse_arguments()
    main(**vars(args))
