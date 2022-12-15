#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Purpose
-------
This module computes the frequency for alleles in a locus
based on allele calling results. It also determines the
number of distinct protein variants and the list of alleles
that code for each protein variant.

Code documentation
------------------
"""


import os
import csv
import argparse

from Bio import SeqIO
from Bio.Seq import Seq


def read_tabular(input_file, delimiter='\t'):
    """ Read TSV file.

    Parameters
    ----------
    input_file : str
        Path to a TSV file.
    delimiter : str
        Delimiter used to separate file columns.

    Returns
    -------
    lines : list
        A list with a sublist per line in the input file.
        Each sublist has the fields that were separated by
        the defined delimiter.
    """

    with open(input_file, 'r') as infile:
        reader = csv.reader(infile, delimiter=delimiter)
        lines = [line for line in reader]

    return lines


def write_to_file(text, output_file, write_mode, end_char):
    """ Writes a single string to a file.

    Parameters
    ----------
    text : str
        A single string to write to the output file.
    output_file : str
        Path to the output file.
    write_mode : str
        Write mode can be 'w', writes text and overwrites
        any text in file, or 'a', appends text to text
        already in file.
    end_char : str
        Character added to the end of the file.
    """

    with open(output_file, write_mode) as out:
        out.write(text+end_char)


def write_lines(lines, output_file):
    """ Writes a list of strings to a file. The strings
        are joined with newlines before being written to
        file.

    Parameters
    ----------
    lines : list
        List with the lines/strings to write to the
        output file.
    output_file : str
        Path to the output file.
    """

    joined_lines = '\n'.join(lines)

    write_to_file(joined_lines, output_file, 'a', '\n')


def translate_sequence(dna_str, table_id):
    """ Translate a DNA sequence using the BioPython package.

    Parameters
    ----------
    dna_str : str
        DNA sequence as string type.
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


def group_by_protein(fasta_file):
    """ Groups DNA sequences based on the protein they code for.

    Parameters
    ----------
    fasta_file : str
        Path to the FASTA file with DNA sequences.

    Returns
    -------
    protein_diversity : dict
        Dictionary with protein sequences as keys and a
        list as value for each key. Each list has the allele
        identifiers of the alleles that code for the protein.
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


def protein_frequencies(protein_groups, locus_classifications):
    """ Creates lines to write Fasta file with proteins
        and determines the frequency of each protein
        based on the alleles that code for that protein.

    Parameters
    ----------
    protein_groups : dict
        Dictionary with protein sequences as keys and
        a list with the IDs of the alleles that code
        for the protein.
    locus_classification : dict
        A dictionary with alleles IDs as keys and a list
        with the sample IDs that have the allele and the
        frequency of the allele based on the results in
        a AlleleCall matrix.

    Returns
    -------
    protein_lines : list
        List with Fasta string records (e.g.: '>1\nMKT').
    prots_freqs : list
        List with one sublist per distinct protein coded
        by the alleles of the locus. Each sublist has
        the protein ID attributed to the protein, the
        frequency of the protein and a list with the IDs
        of the alleles that code for that protein.
    """

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

    return [protein_lines, prots_freqs]


def main(schema, locus_id, allelecall_matrix, output_dir):

    if os.path.isdir(output_dir) is False:
        os.mkdir(output_dir)

    # import AlleleCall matrix
    profiles = read_tabular(allelecall_matrix)

    # get locus column
    try:
        locus_index = profiles[0].index(locus_id+'.fasta')
    except:
    # locus identifier does not have ".fasta" extension
        locus_index = profiles[0].index(locus_id)

    # get locus classifications for all genomes
    locus_classifications = {}
    for line in profiles[1:]:
        allele_id = line[locus_index]
        if 'INF-' in allele_id:
            allele_id = allele_id.replace('INF-', '')

        # get allele identifier and identifiers for samples with that allele
        try:
            allele_id = int(allele_id)
            # do not store missing data that has been converted to "0"
            if allele_id != 0:
                sample_id = line[0]
                locus_classifications.setdefault(allele_id, []).append(sample_id)
        # missing data such as "LNF" cannot be converted to integer
        except ValueError:
            pass

    # determine frequency of each allele
    for k, v in locus_classifications.items():
        v.append(len(v))

    # determine alleles that were not in any genome
    locus_file = os.path.join(schema, locus_id+'.fasta')
    alleles_ids = [int((rec.id).split('_')[-1])
                   for rec in SeqIO.parse(locus_file, 'fasta')]

    # add zero entries for genomes where the locus was not detected
    for a in alleles_ids:
        if a not in locus_classifications:
            locus_classifications[a] = ['', 0]

    # determine frequency at protein level
    protein_groups = group_by_protein(locus_file)

    protein_lines, prots_freqs = protein_frequencies(protein_groups, locus_classifications)

    # write Fasta with proteins
    protein_file = os.path.join(output_dir,
                                '{0}_protein.fasta'.format(locus_id))
    write_lines(protein_lines, protein_file)

    # create output file with protein frequencies
    prots_freqs = sorted(prots_freqs,
                         key=lambda x: x[1],
                         reverse=True)
    prots_freqs_lines = ['{0}\t{1}\t{2}'.format(l[0], l[1], ','.join(l[2]))
                         for l in prots_freqs]
    prots_freqs_header = 'protein_id\tfrequency\talleles_ids'
    prots_freqs_lines = [prots_freqs_header] + prots_freqs_lines
    protein_freqs_file = os.path.join(output_dir,
                                      '{0}_protein_frequencies.tsv'.format(locus_id))
    write_lines(prots_freqs_lines, protein_freqs_file)

    # create output file with DNA sequences frequencies
    output_lines = [[k, v[-1], v[:-1]]
                    for k, v in locus_classifications.items()]
    # sort based on decreased frequency
    output_lines = sorted(output_lines,
                          key=lambda x: x[1],
                          reverse=True)
    output_header = 'allele_id\tfrequency\tsamples'
    output_text = ['{0}\t{1}\t{2}'.format(l[0], l[1], ','.join(l[2]))
                   for l in output_lines]
    output_text = [output_header] + output_text

    output_file = os.path.join(output_dir,
                               '{0}_frequencies.tsv'.format(locus_id))
    write_lines(output_text, output_file)


def parse_arguments():

    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument('-s', type=str, required=True,
                        dest='schema',
                        help='Path to the schema directory.')

    parser.add_argument('-l', type=str, required=True,
                        dest='locus_id',
                        help='Locus identifier (without ".fasta" extension).')

    parser.add_argument('-m', type=str, required=True,
                        dest='allelecall_matrix',
                        help='Path to file that contains a matrix with '
                             'allele calling results.')

    parser.add_argument('-o', type=str, required=True,
                        dest='output_dir',
                        help='Path to output directory.')

    args = parser.parse_args()

    return args


if __name__ == '__main__':

    args = parse_arguments()

    main(**vars(args))
