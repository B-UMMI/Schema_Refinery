#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Purpose
-------

This sub-module aligns schema representative sequences
against records in Genbank files to extract relevant
annotations.

Code documentation
------------------
"""

import os
import csv

from Bio import SeqIO

try:
    from utils import blast_functions as bf
    from utils.sequence_functions import translate_sequence
    from utils import file_functions as ff
except ModuleNotFoundError:
    from SchemaRefinery.utils import blast_functions as bf
    from SchemaRefinery.utils.sequence_functions import translate_sequence
    from SchemaRefinery.utils import file_functions as ff


def get_protein_annotation_fasta(seqRecord):
    """Get the translated protein from a Genbank file.

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

    Notes
    -----
    Source: https://github.com/LeeBergstrand/Genbank-Downloaders/blob/d904c92788696b02d9521802ebf1fc999a600e1b/SeqExtract.py#L48
    """
    fasta = []
    fasta_dict = {}
    features = seqRecord.features  # Each sequence has a list (called features) that stores seqFeature objects.
    for feature in features:  # For each feature on the sequence
        if feature.type == "CDS":  # CDS means coding sequence (These are the only features we're interested in)
            featQualifiers = feature.qualifiers  # Each feature contains a dictionary called qualifiers which contains
            # data about the sequence feature (for example the translation)

            # Gets the required qualifier. Uses featQualifers.get to return the qualifier or a default value if the quantifier
            # is not found. Calls strip to remove unwanted brackets and ' from qualifier before storing it as a string.
            protein_id = str(featQualifiers.get('protein_id', '')).strip('\'[]')

            if protein_id == 'no_protein_id':
                continue  # Skips the iteration if protein has no id.

            gene = str(featQualifiers.get('gene', '')).strip('\'[]')
            product = str(featQualifiers.get('product', 'no_product_name')).strip('\'[]')
            translated_protein = str(featQualifiers.get('translation', 'no_translation')).strip('\'[]')

            fasta.append(('>' + protein_id + '|' + gene + '|' + product + '\n' + translated_protein))
            fasta_dict[translated_protein] = [protein_id, product, gene]

    return fasta, fasta_dict


def genbank_annotations(genbank_files: str, schema_directory: str,
                        output_directory: str, cpu_cores: int,
                        bsr: float):

    output_directory = os.path.join(output_directory, 'genbank_annotations')
    ff.create_directory(output_directory)

    gbk_files = [os.path.join(genbank_files, f)
                 for f in os.listdir(genbank_files)]
    gbk_files.sort()

    fasta = []
    fasta_dict = {}
    for f in gbk_files:
        recs = [rec for rec in SeqIO.parse(f, 'genbank')]
        for r in recs:
            outl, outd = get_protein_annotation_fasta(r)
            fasta.extend(outl)
            fasta_dict.update(outd)

    selected_file = os.path.join(output_directory, 'selected_cds.fasta')
    with open(selected_file, 'w') as outfile:
        fasta_text = '\n'.join(fasta)
        outfile.write(fasta_text)

    # BLAST alleles for each locus against file with all CDSs from origin genomes
    reps_dir = os.path.join(schema_directory, 'short')
    rep_files = [os.path.join(reps_dir, f)
                 for f in os.listdir(reps_dir)
                 if f.endswith('.fasta')]

    # Get all representative sequences into same file
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

    # Save new reps into same file
    prot_file = os.path.join(output_directory, 'reps_prots.fasta')
    with open(prot_file, 'w') as pinfile:
        pinfile.write('\n'.join(reps))

    # Create BLASTdb
    blastdb_path = os.path.join(output_directory, 'reps_db')
    bf.make_blast_db(prot_file, blastdb_path, 'prot')

    blastout = os.path.join(output_directory, 'blastout.tsv')
    bf.run_blast('blastp', blastdb_path, selected_file, blastout,
                 max_hsps=1, threads=cpu_cores, max_targets=1)

    # Import BLAST results
    with open(blastout, 'r') as at:
        blast_results = list(csv.reader(at, delimiter='\t'))

    # Convert ids
    for r in blast_results:
        r[1] = reps_ids[int(r[1])]

    # Create mapping between sequence IDs and sequence
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

    # Get identifiers mapping
    ids_to_name = {v[0]: v[1:] for k, v in fasta_dict.items()}

    # Add names
    for k in best_matches:
        best_matches[k].extend(ids_to_name[best_matches[k][0]])

    # Concatenate reps and get self-score
    reps_blast_out = os.path.join(output_directory, 'concat_reps_self.tsv')
    bf.run_blast('blastp', blastdb_path, prot_file, reps_blast_out,
                 max_hsps=1, threads=cpu_cores, ids_file=None, blast_task=None,
                 max_targets=1)

    # Import self results
    with open(reps_blast_out, 'r') as at:
        reps_blast_results = list(csv.reader(at, delimiter='\t'))

    # Convert ids
    for r in reps_blast_results:
        r[0] = reps_ids[int(r[0])]
        r[1] = reps_ids[int(r[1])]

    reps_scores = {l[0]: l[-1] for l in reps_blast_results}

    for k in best_matches:
        best_matches[k].append(float(best_matches[k][1])/float(reps_scores[k]))

    final_best_matches = {k: v for k, v in best_matches.items() if v[5] >= bsr}
    print('Extracted annotations for {0} loci.'.format(len(final_best_matches)))

    # Save annotations
    header = 'Locus\tgenebank_origin_id\tgenebank_origin_product\tgenebank_origin_name\tgenebank_origin_bsr'
    annotations_file = os.path.join(output_directory, 'genbank_annotations.tsv')
    with open(annotations_file, 'w') as at:
        outlines = [header] + ['{0}\t{1}\t{2}\t{3}\t{4}'.format(k.split("_")[0], v[0], v[3], v[4], v[5]) for k, v in final_best_matches.items()]
        outtext = '\n'.join(outlines)
        at.write(outtext+'\n')

    return annotations_file
