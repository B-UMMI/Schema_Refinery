#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Purpose
-------
This script determines the start and stop positions in the
genome of origin for the first representative sequences in
a schema.

Code documentation
------------------
"""


import os
import argparse

from Bio import SeqIO, Seq


def main(input_files, schema_dir, output_dir):

    if os.path.isdir(output_dir) is False:
        os.mkdir(output_dir)

    # dict with genome assembly identifier to full path
    genomes = {os.path.basename(f): os.path.join(input_files, f)
               for f in os.listdir(input_files)}
    # replace characters to get same identiifers as in schema
    genomes = {k.replace('_', '-'): v for k, v in genomes.items()}
    genomes = {k.split('.fasta')[0]: v for k, v in genomes.items()}

    # dict with locus identifier to full path
    genes = {f.split('.fasta')[0]: os.path.join(schema_dir, f)
             for f in os.listdir(schema_dir)
             if '.fasta' in f}

    # dict with locus identifier to representative sequence
    genes_seqs = {}
    for gene, gene_path in genes.items():
        # only get first sequence/representative
        rep_record = SeqIO.parse(gene_path, 'fasta').__next__()
        rep_seq = (str(rep_record.seq), rep_record.id)
        genes_seqs[gene] = rep_seq

    # dict with locus identifier to genome of origin
    genes_genomes = {}
    for gene, gene_path in genes.items():
        gene_genome = gene.split('-protein')[0]
        for genome, genome_path in genomes.items():
            if gene_genome in genome:
                genes_genomes[gene] = genome_path

    # group genes by genome
    genes_groups = {}
    for k, v in genes_genomes.items():
        genes_groups.setdefault(v, []).append(k)

    # get start and stop positions
    missing = []
    genes_positions = {}
    for genome, gene_group in genes_groups.items():
        genome_contigs = {str(rec.seq): rec.id
                          for rec in SeqIO.parse(genome, 'fasta')}
        for g in gene_group:
            seq = genes_seqs[g][0]
            for c in genome_contigs:
                sense_strand = c
                contig_id = genome_contigs[c]
                if seq in sense_strand:
                    start = sense_strand.index(seq)
                    genes_positions[g] = [start, start+len(seq), 1, contig_id]
                else:
                    rev_seq = Seq.reverse_complement(seq)
                    if rev_seq in sense_strand:
                        start = sense_strand.index(rev_seq)
                        genes_positions[g] = [start, start+len(rev_seq), 0, contig_id]

            if g not in genes_positions:
                missing.append(g)

    print('Could not determine start and stop positions for '
          '{0} loci.'.format(missing))

    # save results
    header = 'locus\tcontig\tstart\tstop\tstrand'
    positions_outfile = os.path.join(output_dir, 'loci_positions.tsv')
    with open(positions_outfile, 'w') as outfile:
        outlines = [header] + ['{0}\t{1}\t{2}\t{3}\t{4}'.format(k, v[-1], v[0], v[1], v[2])
                               for k, v in genes_positions.items()]
        outtext = '\n'.join(outlines)
        outfile.write(outtext)


def parse_arguments():

    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument('-i', type=str, required=True,
                        dest='input_files',
                        help='Path to the directory that contains the '
                             'genome assemblies in Fasta format.')

    parser.add_argument('-s', type=str, required=True,
                        dest='schema_dir',
                        help='Path to the schemas\'s directory.')

    parser.add_argument('-o', type=str, required=True,
                        dest='output_dir',
                        help='Path to the output directory.')

    args = parser.parse_args()

    return args


if __name__ == '__main__':

    args = parse_arguments()
    main(**vars(args))
