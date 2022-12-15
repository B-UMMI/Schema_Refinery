#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Purpose
-------
This script aligns schema representative sequences
against records in Genbank files to extract relevant
annotations.

Code documentation
------------------
"""

import os
import csv
import shutil
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
    """ Determines the reverse complement of given DNA strand.

    Parameters
    ----------
    dna_sequence : str
        String representing a DNA sequence.

    Returns
    -------
    reverse_complement_strand : str
        The reverse complement of the DNA sequence,
        without lowercase letters.
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


def getProteinAnnotationFasta(seqRecord):
    """ Gets the translated protein from a Genbank file.
    Source: https://github.com/LeeBergstrand/Genbank-Downloaders/blob/d904c92788696b02d9521802ebf1fc999a600e1b/SeqExtract.py#L48
    
    Parameters
    ----------
    seqRecord : Biopython SeqRecord
        BioPython sequence record object.

    Returns
    -------
    fasta : list
        List containing the protein in fasta format.
    fasta_dict : dict
        Dict containing the translated protein as key and the values are
        list containing the protein_id, the product and the gene name.

    """
    fasta = []
    fasta_dict = {}
    features = seqRecord.features  # Each sequence has a list (called features) that stores seqFeature objects.
    for feature in features:  # For each feature on the sequence
        
        if feature.type == "CDS":  # CDS means coding sequence (These are the only feature we're interested in)
            featQualifers = feature.qualifiers  # Each feature contains a dictionary called qualifiers which contains
            # data about the sequence feature (for example the translation)

            # Gets the required qualifier. Uses featQualifers.get to return the qualifier or a default value if the quantifier
            # is not found. Calls strip to remove unwanted brackets and ' from qualifier before storing it as a string.
            protein_id = str(featQualifers.get('protein_id', 'no_protein_id')).strip('\'[]')
            
            if protein_id == 'no_protein_id':
                continue  # Skips the iteration if protein has no id.
            
            gene = str(featQualifers.get('gene', '')).strip('\'[]')
            product = str(featQualifers.get('product', 'no_product_name')).strip('\'[]')
            translated_protein = str(featQualifers.get('translation', 'no_translation')).strip('\'[]')
            
            fasta.append((">" + protein_id + "|" + gene + "|" + product + "\n" + translated_protein))
            fasta_dict[translated_protein] = [protein_id, product, gene]
            
    return fasta, fasta_dict


def main(input_files, schema_dir, output_dir, cpu_cores):

    if os.path.isdir(output_dir) is False:
        os.mkdir(output_dir)

    gbk_files = [os.path.join(input_files, f) for f in os.listdir(input_files)]
    
    gbk_files.sort()
        
    fasta = []
    fasta_dict = {}
    for f in gbk_files:
        recs = [rec for rec in SeqIO.parse(f, 'genbank')]
        
        for r in recs:
            outl, outd = getProteinAnnotationFasta(r)
            fasta.extend(outl)
            fasta_dict.update(outd)

    selected_file = os.path.join(output_dir, 'selected_cds.fasta')
    with open(selected_file, 'w') as outfile:
        fasta_text = '\n'.join(fasta)
        outfile.write(fasta_text)
            

    # BLAST alleles for each locus against file with all CDSs from origin genomes
    reps_dir = os.path.join(schema_dir, 'short')
    rep_files = [os.path.join(reps_dir, f) for f in os.listdir(reps_dir) if f.endswith('.fasta') is True]

    # get all representative sequences into same file
    reps = []
    reps_ids = {}
    start = 1
    for f in rep_files:
        first_rep =  SeqIO.parse(f, 'fasta').__next__()
        prot = translate_sequence(str(first_rep.seq), 11)
        sequence = '>{0}\n{1}'.format(start, prot)
        reps_ids[start] = first_rep.id
        reps.append(sequence)
        start += 1

    # save new reps into same file
    prot_file = os.path.join(output_dir, 'reps_prots.fasta')
    with open(prot_file, 'w') as pinfile:
        pinfile.write('\n'.join(reps))

    # create BLASTdb
    blastdb_path = os.path.join(output_dir, 'reps_db')
    make_blast_db('makeblastdb', prot_file, blastdb_path, 'prot')

    blastout = os.path.join(output_dir, 'blastout.tsv')
    run_blast('blastp', blastdb_path, selected_file, blastout,
              max_hsps=1, threads=cpu_cores, ids_file=None,
              blast_task=None, max_targets=1, ignore=None)

    # import BLAST results
    with open(blastout, 'r') as at:
        blast_results = list(csv.reader(at, delimiter='\t'))

    # reverse ids
    for r in blast_results:
        r[1] = reps_ids[int(r[1])]

    # create mapping between sequence IDs and sequence
    selected_inverse = {v[0]: [k, v[2]] for k, v in fasta_dict.items()}

    best_matches = {}
    for rec in blast_results:
        try:
            query = rec[0].split('|')[0]
            subject = rec[1]
            score = rec[-1]
            query_name = selected_inverse[query][1]
        except KeyError:
            continue
        if subject in best_matches:
            if best_matches[subject][2] == '' and query_name != 'NA':
                best_matches[subject] = [query, score, query_name]
            elif best_matches[subject][2] != '' and query_name != 'NA':
                if float(score) > float(best_matches[subject][1]):
                    best_matches[subject] = [query, score, query_name]
            elif best_matches[subject][2] == '' and query_name == 'NA':
                if float(score) > float(best_matches[subject][1]):
                    best_matches[subject] = [query, score, query_name]
        else:
            best_matches[subject] = [query, score, query_name]

    # get identifiers mapping
    ids_to_name = {v[0]: v[1:] for k, v in fasta_dict.items()}
    
    # add names
    for k in best_matches:
        best_matches[k].extend(ids_to_name[best_matches[k][0]])

    # concatenate reps and get self-score
    reps_blast_out = os.path.join(output_dir, 'concat_reps_self.tsv')
    run_blast('blastp', blastdb_path, prot_file, reps_blast_out,
              max_hsps=1, threads=cpu_cores, ids_file=None, blast_task=None,
              max_targets=1, ignore=None)

    # import self_results
    with open(reps_blast_out, 'r') as at:
        reps_blast_results = list(csv.reader(at, delimiter='\t'))

    # reverse ids
    for r in reps_blast_results:
        r[0] = reps_ids[int(r[0])]
        r[1] = reps_ids[int(r[1])]

    reps_scores = {l[0]: l[-1] for l in reps_blast_results}

    for k in best_matches:
        best_matches[k].append(float(best_matches[k][1])/float(reps_scores[k]))

    final_best_matches = {k: v for k, v in best_matches.items() if v[5] >= 0.6}
    print('Extracted annotations for {0} loci.'.format(len(final_best_matches)))

    # save annotations
    header = 'Locus_ID\torigin_id\torigin_product\torigin_name\torigin_bsr'
    annotations_file = os.path.join(output_dir, 'genbank_annotations.tsv')
    with open(annotations_file, 'w') as at:
        outlines = [header] + ['{0}\t{1}\t{2}\t{3}\t{4}'.format(k.split("_")[0], v[0], v[3], v[4], v[5]) for k, v in final_best_matches.items()]
        outtext = '\n'.join(outlines)
        at.write(outtext+'\n')

    print('Results available at {0}'.format(annotations_file))


def parse_arguments():

    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument('-i', type=str, required=True,
                        dest='input_files',
                        help='Path to the directory that contains '
                             'Genbank files with annotations to '
                             'extract.')

    parser.add_argument('-s', type=str, required=True,
                        dest='schema_dir',
                        help='Path to the schema\'s directory.')

    parser.add_argument('-o', type=str, required=True,
                        dest='output_dir',
                        help='Path to the output directory.')

    parser.add_argument('-cpu', type=int, required=False,
                        default=1, dest='cpu_cores',
                        help='Number of CPU cores used by BLAST.')

    args = parser.parse_args()

    return args


if __name__ == '__main__':

    args = parse_arguments()
    main(**vars(args))
