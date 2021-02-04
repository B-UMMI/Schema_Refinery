#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

@author: rfm
"""


import os
import csv
import argparse

from Bio import SeqIO
from Bio.Seq import Seq


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


def group_by_protein(fasta_file):
    """ Groups DNA sequences based on the protein they code for.

        Parameters
        ----------
        fasta_file : str
            Path to the FASTA file with DNA sequences.

        Returns
        -------
        protein_diversity : dict
            Dictionary with a gene identifier as key and
            another dictionary as value. The nested dictionary
            has protein sequences as keys and a list as value
            for each key. Each list has the allele identifiers
            and sequences that code for that protein, organized
            in tuples.
    """

    protein_diversity = {}
    for record in SeqIO.parse(fasta_file, 'fasta'):
        seqid = record.id
        allele_id = seqid.split('_')[-1]
        sequence = str(record.seq)
        try:
            protein = str(translate_sequence(sequence, 11))
        except Exception:
            continue

        protein_diversity.setdefault(protein, []).append(allele_id)

    return protein_diversity


def main(schema, locus_id, allelecall_matrix, output_dir):

    if os.path.isdir(output_dir) is False:
        os.mkdir(output_dir)

    # import AlleleCall matrix
    with open(allelecall_matrix, 'r') as infile:
        profiles = list(csv.reader(infile, delimiter='\t'))

    # get locus column
    locus_index = profiles[0].index(locus_id+'.fasta')

    # get locus classifications for all genomes
    locus_classifications = {}
    for line in profiles[1:]:
        allele_id = line[locus_index]
        if 'INF-' in allele_id:
            allele_id = allele_id.replace('INF-', '')

        try:
            allele_id = int(allele_id)
            sample_id = '.'.join(line[0].split('.')[:-1])
            locus_classifications.setdefault(allele_id, []).append(sample_id)
        except ValueError:
            pass

    # determine frequency of each allele
    for k, v in locus_classifications.items():
        v.append(len(v))

    # determine alleles that were not in any genome
    locus_file = os.path.join(schema, locus_id+'.fasta')
    alleles_ids = [int((rec.id).split('_')[-1])
                   for rec in SeqIO.parse(locus_file, 'fasta')]

    for a in alleles_ids:
        if a not in locus_classifications:
            locus_classifications[a] = ['', 0]

    # determine frequency at protein level
    protein_groups = group_by_protein(locus_file)

    # write file with proteins and add protein freqs
    protid = 1
    protein_lines = []
    prots_freqs = []
    for k, v in protein_groups.items():
        line = '>{0}\n{1}'.format(protid, k)
        protein_lines.append(line)
        total_freq = 0
        for a in v:
            total_freq += locus_classifications[int(a)][-1]
        prots_freqs.append([protid, total_freq, v])
        protid += 1

    # write Fasta with proteins
    protein_file = os.path.join(output_dir, '{0}_protein.fasta'.format(locus_id))
    with open(protein_file, 'w') as outfile:
        outfile.write('\n'.join(protein_lines))

    # create output file with protein frequencies
    prots_freqs = sorted(prots_freqs, key=lambda x: x[1], reverse=True)
    prots_freqs_lines = ['{0}\t{1}\t{2}'.format(l[0], l[1], ','.join(l[2]))
                         for l in prots_freqs]
    prots_freqs_header = 'protein_id\tfrequency\talleles_ids'
    prots_freqs_lines = [prots_freqs_header] + prots_freqs_lines
    protein_freqs_file = os.path.join(output_dir,
                                      '{0}_protein_frequencies.tsv'.format(locus_id))
    with open(protein_freqs_file, 'w') as outfile:
        outfile.write('\n'.join(prots_freqs_lines))

    # create output file with DNA sequences frequencies
    output_lines = [[k, v[-1], v[:-1]] for k, v in locus_classifications.items()]
    # sort based on decreased frequency
    output_lines = sorted(output_lines, key=lambda x: x[1], reverse=True)
    output_header = 'allele_id\tfrequency\tsamples'
    output_text = ['{0}\t{1}\t{2}'.format(l[0], l[1], ','.join(l[2]))
                   for l in output_lines]
    output_text = [output_header] + output_text
    output_text = '\n'.join(output_text)

    output_file = os.path.join(output_dir, '{0}_frequencies.tsv'.format(locus_id))
    with open(output_file, 'w') as outfile:
        outfile.write(output_text)


def parse_arguments():

    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument('-s', type=str, required=True,
                        dest='schema',
                        help='')

    parser.add_argument('-l', type=str, required=True,
                        dest='locus_id',
                        help='')

    parser.add_argument('-m', type=str, required=True,
                        dest='allelecall_matrix',
                        help='')

    parser.add_argument('-o', type=str, required=True,
                        dest='output_dir',
                        help='')

    args = parser.parse_args()

    return args


if __name__ == '__main__':

    args = parse_arguments()

    main(**vars(args))
