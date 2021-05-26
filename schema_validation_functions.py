#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Purpose
-------

This module has the collection of functions that are used to validate
schemas. The schema structure should respect a set of criteria that are
used during the CreateSchema and AlleleCall processes. All alleles for
a gene present in the schema should be selected only if they share a
BSR >= 0.6 with at least one representative allele of that gene.
Representative alleles from different loci should not share a BSR >= 0.6.
However, the process does not guarantee that all alleles for a gene in
the schema share a BSR >= 0.6 or that non-representative sequences from
one gene cannot have a BSR >= 0.6 with alleles from other gene. The
functions in this module were created with the intent of validating
schemas and help determine if a given schema or subset of a schema
respects the necessary criteria and does not deviate too much from
those criteria due to aspects that are not controlled during the
CreateSchema or AlleleCall processes.

Several functions in this module require that 'Biopython' be installed
in the working environment.
"""


import os
import csv
import time
import shutil
import itertools
import traceback
import subprocess
from collections import Counter
from multiprocessing import Pool

from Bio import SeqIO
from Bio.Seq import Seq


def schema_dictionary(schema_files, output_directory):
    """
    """

    loci_reassigned = {}
    # for each FASTA file that has the alleles for a locus
    for locus, file in schema_files.items():
        records = [((rec.id).split('_')[-1], str(rec.seq))
                   for rec in SeqIO.parse(file, 'fasta')]
        fasta_records = ['>{0}_{1}\n{2}'.format(locus, rec[0], rec[1])
                         for rec in records]
        fasta_text = '\n'.join(fasta_records)
        output_file = os.path.join(output_directory, os.path.basename(file))
        loci_reassigned[locus] = output_file

        with open(output_file, 'w') as outfile:
            outfile.write(fasta_text+'\n')

    return loci_reassigned


def translate_schema(fasta_files, output_directory):
    """ Creates a dictionary that stores the protein sequences
        coded in the DNA sequences of each allele per gene in
        the schema.

        Parameters
        ----------
        fasta_files : str
            

        Returns
        -------
        protein_schema_dict : dict
            A dictionary with gene identifiers as keys. Each key
            has a sub dictionary with alleles identifiers as keys
            and alleles protein sequences as values.
    """

    protein_files = {}
    for locus, file in fasta_files.items():
        records = [(rec.id, str(rec.seq))
                   for rec in SeqIO.parse(file, 'fasta')]
        protein_records = [(rec[0], str(translate_sequence(rec[1], 11)))
                           for rec in records]

        protein_fasta = ['>{0}\n{1}'.format(rec[0], rec[1])
                         for rec in protein_records]

        protein_text = '\n'.join(protein_fasta)
        output_file = os.path.join(output_directory,
                                   os.path.basename(file)+'_protein')
        protein_files[locus] = output_file
        with open(output_file, 'w') as outfile:
            outfile.write(protein_text+'\n')

    return protein_files


def translate_sequence(dna_str, table_id):
    """ Translates a DNA sequence using the BioPython package.

        Parameters
        ----------
        dna_str : str
            String representing a DNA sequence.
        table_id : int
            Translation table identifier.

        Returns
        -------
        protseq : Bio.Seq.Seq
            Protein sequence created by translating the
            input DNA sequence.
    """

    myseq_obj = Seq(dna_str)
    protseq = Seq.translate(myseq_obj, table=table_id, cds=True)

    return protseq


def decode_str(str_list, encoding):
    """ Decodes bytes objects in the input list and
        strips decoded strings from whitespaces and
        newlines.

        Parameters
        ----------
        str_list
            List with string or bytes objects to decode
            and strip of whitespaces and newlines.
        encoding : str
            Encoding codec to use.

        Returns
        -------
        decoded : list
            List with strings without whitespaces or
            newlines.
    """

    decoded = [m.decode(encoding).strip()
               if type(m) == bytes
               else m.strip()
               for m in str_list]

    return decoded


def filter_list(lst, remove):
    """ Removes elements from a list.

        Parameters
        ----------
        lst : list
            Input list.
        remove : list
            List of elements to remove from input list.

        Returns
        -------
        filtered_list : list
            List without the removed elements.
    """

    filtered_list = list(set(lst) - set(remove))

    return filtered_list


def make_blast_db(makeblastdb_path, input_fasta, output_path, db_type,
                  ignore=None):
    """ Creates a BLAST database.

        Parameters
        ----------
        makeblastdb_path : str
            Path to the 'makeblastdb' executable.
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
        ignore : list of None
            List with BLAST warnings that should be ignored.

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

    if len(stderr) > 0:
        stderr = decode_str(stderr, 'utf8')
        if ignore is not None:
            stderr = filter_list(stderr, ignore)

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
            Maximum number of target/subject sequences
            to align against.
        ignore : list or None
            List with BLAST warnings that should be ignored.

        Returns
        -------
        stderr : str
            String with errors raised during BLAST execution.
    """

    blast_args = [blast_path, '-db', blast_db, '-query', fasta_file,
                  '-out', blast_output, '-outfmt', '6 qseqid sseqid score',
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

    if len(stderr) > 0:
        stderr = decode_str(stderr, 'utf8')
        if ignore is not None:
            stderr = filter_list(stderr, ignore)

    return stderr


def read_blast_tabular(blast_tabular_file):
    """ Read a file with BLAST results in tabular format

        Args: 
            blast_tabular_file (str): path to output file of BLAST.

        Returns:
            blasting_results (list): a list with a sublist per line in the input
            file.
    """

    with open(blast_tabular_file, 'r') as blastout:
        blasting_results = []
        reader = csv.reader(blastout, delimiter='\t')
        for row in reader:
            blasting_results.append(row)

    return blasting_results


def create_blast_results_dictionary(blast_results):
    """ Creates a dictionary based on the BLAST results of matches given as input.

        Args:
            blast_results (list): list with a sublist for each line that was in the 
            BLAST output file.

        Returns:
            blast_matches (dict): a dictionary with genes/loci identifiers as keys.
            Values are dictionaries with alleles identifiers as keys and dictionaries with 
            hits identifiers as keys and a list of BLASTp information about the hits as values.
    """

    blast_matches = {}
    # for every hit in the BLAST results
    for match in blast_results:

        # get unique name for the query and hit
        query_name = match[0].strip()
        hit_name = match[1].strip()

        gene_id = query_name.split('_')[0]

        # get the relevant information
        identity = match[2]
        query_length = match[3]
        query_start = match[4]
        query_end = match[5]
        subject_start = match[6]
        subject_end = match[7]
        #evalue = match[8]
        bitscore = match[9]
        raw_score = match[10]

        # create the gene/locus key with the first query and hit information
        if gene_id not in blast_matches:
            blast_matches[gene_id] = {query_name:{hit_name: [identity, query_length, query_start,
                                                                  query_end, subject_start, subject_end, 
                                                                  bitscore, raw_score]}}
        
        # add query and hit information if the gene/locus is in the dictionary
        elif gene_id in blast_matches:
            # if the query already has hits, just add the hit
            if query_name in blast_matches[gene_id]:
                blast_matches[gene_id][query_name][hit_name] = [identity, query_length, query_start, 
                                                                     query_end, subject_start, subject_end, 
                                                                     bitscore, raw_score]

            # add the query key and the information about the first hit
            elif query_name not in blast_matches[gene_id]:
                blast_matches[gene_id][query_name] = {hit_name: [identity, query_length, query_start, 
                                                                      query_end, subject_start, subject_end, 
                                                                      bitscore, raw_score]}
    
    return blast_matches


def complete_loci_matches_info(blast_matches, protein_schema, representatives):
    """ Creates a new dictionary with more information about the BLAST results.
    
        Args:
            blast_matches (dict): dictionary with the information about the hits for the 
            alleles of one gene/locus. Only includes the fields present in the BLAST output.
            protein_schema (dict): dictionary with genes/loci identifiers as keys and 
            protein sequences as values.
            representatives (list): list with the identifiers of the representative alleles 
            of all genes/loci.
        
        Returns:
            matches_processed_info (dict): a dictionary similar to 'blast_matches' but with 
            additional information about the BLAST matches.
        
        Example:
            >>> complete_loci_matches_info()
    """

    # determine cases of different gene alleles with BSR >= 0.6
    matches_complete_info = {}
    distance_matrix = {}
    # for each gene prefix
    for gene in blast_matches:
        # get the representative alleles of the gene
        alleles = blast_matches[gene]
        matches_complete_info[gene] = {}
        distance_matrix[gene] = {}

        # get query information and calculate length mode and ASM and ALM thresholds
        for allele in alleles:
            current_allele = alleles[allele]
            
            self_blast_score = current_allele[allele][7]
            distance_matrix[gene][allele] = {}
            
            matches_complete_info[gene][allele] = {}

            query_proteins = protein_schema[int(gene)]
            query_info = [[query_proteins[p],len(query_proteins[p])] for p in query_proteins if allele in p][0]

            # determine the locus length mode
            locus_length_mode = calculate_locus_mode(query_proteins)
            asm_factor, alm_factor = relative_to_mode(locus_length_mode, 0.2)
            
            # for each hit found for the query
            for hit, stats in current_allele.items():

                blast_score_ratio = calculate_blast_score_ratio(stats[7], self_blast_score)

                distance_matrix[gene][allele][hit] = blast_score_ratio
                
                # get hit protein sequences
                hit_gene = int(hit.split('_')[0])
                hit_proteins = protein_schema[hit_gene]

                # determine hit protein length
                hit_info = [[hit_proteins[p],len(hit_proteins[p])] for p in hit_proteins if hit in p][0]
                    
                # determine if it might be an ASM or ALM case
                mode_deviation = determine_mode_deviation(hit_info[1], asm_factor, alm_factor)
            
                # determine if query is contained in hit
                # should only happen in cases with high BSR
                contained = determine_contained(query_info[0], hit_info[0])
                
                match_type = '{0} & {1}'.format(mode_deviation, contained)
                
                representativeness_string = match_representatives(allele, hit, representatives)
                
                # add the match with all information to the dictionary
                matches_complete_info[gene][allele][hit] = [allele, hit] + stats + [str(self_blast_score), 
                                                              str(blast_score_ratio), str(query_info[1]), str(hit_info[1]), 
                                                              str(locus_length_mode), match_type, representativeness_string]

    return [matches_complete_info, distance_matrix]


def match_representatives(query_id, hit_id, representatives):
    """ Identifies sequences that are representative alleles.
    
        Args:
            query_id (str): identifier of the query sequence.
            hit_id (str): identifier of the subject sequence.
            representatives (list): list with all identifiers of representative
            alleles.
            
        Returns:
            query_rep + hit_rep (str): a string that has two characters, one for 
            the query sequence and another one for the subject sequence. '1' indicates 
            that the sequence is a representative allele, '0' otherwise.
    """

    query_rep = '1' if query_id in representatives else '0'
        
    hit_rep = '1' if hit_id in representatives else '0'
    
    return query_rep + hit_rep


def determine_contained(query_protein, subject_protein):
    """ Determines if the query or the subject are equal or if any of those is 
        contained in the other.
        
        Args:
            query_protein (str): protein sequence of the query.
            subject_protein (str): protein sequence of the subject.
        
        Returns:
            contained (str): string that indicates if the sequences are equal or if 
            any of them contains the other.
    """
    
    if query_protein == subject_protein:
        contained = 'EQUAL'
    elif query_protein in subject_protein:
        contained = 'QUERY IN SUBJECT'
    elif subject_protein in query_protein:
        contained = 'SUBJECT IN QUERY'
    else:
        contained = 'NOT CONTAINED'

    return contained


def determine_mode_deviation(protein_length, asm_factor, alm_factor):
    """ Determines if a given protein is too short or too long according to its 
        length and the matching allele length mode.
        
        Args:
            protein_length (int): length of the given protein.
            asm_factor (int): protein length value that indicates that a protein is 
            too short.
            alm_factor (int): protein length value that indicates that a protein is 
            too long.
        
        Returns:
            mode_deviation (str): string indicating if the given protein length 
            is in the mode value interval that is accepted or if the protein length 
            is an ASM or ALM.
    """
    
    if protein_length <= asm_factor:
        mode_deviation = 'ASM'
    elif protein_length >= alm_factor:
        mode_deviation = 'ALM'
    else:
        mode_deviation = 'IN MODE'

    return mode_deviation


def calculate_blast_score_ratio(hit_score, self_score):
    """ Calculates the BLAST Score Ratio for a given match.
        
        Args:
            hit_score (float): the BLAST raw score for the match between the query 
            protein and the subject protein.
            self_score (float): the BLAST raw score of the match between the query 
            and itself.
        
        Returns:
            blast_score_ratio (float): the BLAST Score Ratio for the given match.
    """
    
    blast_score_ratio = float(hit_score) / float(self_score)
    
    return blast_score_ratio


def calculate_locus_mode(protein_dictionary):
    """ Calculates the locus length mode.
    
        Args:
            protein_dictionary (dict): dictionary with sequence identifiers as keys 
            and protein sequences as values.
                
        Returns:
            locus_length_mode (int): the mode of the protein sequence length of the 
            given gene/locus.
    """
    
    allele_sizes = []
    for p in protein_dictionary:
        protein_length = len(protein_dictionary[p])
        allele_sizes.append(protein_length)
    
    # return most common length value
    locus_length_mode = Counter(allele_sizes).most_common()[0][0]
    
    return locus_length_mode


def relative_to_mode(locus_length_mode, size_threshold=0.2):
    """ Determine the sequence length threshold that will be used to define 
        ASM and ALM cases.
        
        Args:
            locus_length_mode (int): sequence length mode of the gene/locus.
            size_threshold (float): the ratio used to define the acceptable threshold 
            for ASM and ALM cases.
    """
    
    allele_smaller_than_mode = locus_length_mode - (locus_length_mode * size_threshold)
    allele_larger_than_mode = locus_length_mode + (locus_length_mode * size_threshold)
    
    return [allele_smaller_than_mode, allele_larger_than_mode]


def create_within_loci_lines(matches, blast_score_ratio, header, genes_ids, alleles_ids):
    """ Creates the lines for all matches between alleles of the same locus and 
        a separate list with just the matches between alleles of the same locus with 
        a BSR < 0.6.
        
        Args:
            matches (list): list with one dictionary per locus.
            blast_score_ratio (float): BLAST Score Ratio value.
            header (str): header for the files with all matches and with just 
            the low BSR matches.
        
        Returns:
            all_matches_lines (list): list where the first element is the header 
            for the file and all remaining lines are the matches found by BLAST, 
            without restrictions considering the BSR value.
            low_bsr_lines (list): list where the first element is the header for 
            the file with the matches that have a BSR value below the considered 
            threshold.
    """
    
    # append header string/line
    all_matches_lines = []
    all_matches_lines.append(header)
    
    low_bsr_lines = []
    low_bsr_lines.append(header)
    
    # for each locus
    for dictionary in matches:
        for locus in dictionary:
            apppended = 0
            gene_id = genes_ids[int(locus)]
            # get the representatives of the locus
            representatives = dictionary[locus]
            for rep in representatives:
                # get the hits for the current representative
                rep_id = alleles_ids[rep]
                hits = representatives[rep]
                # join match information into single string
                for h in hits:
                    hits[h][0] = rep_id
                    hit_id = alleles_ids[h]
                    hits[h][1] = hit_id
                    new_line = ','.join(hits[h])
                    
                    all_matches_lines.append(new_line+'\n')
                    
                    # add string to the list of matches with BSR lower than the threshold
                    if float(hits[h][11]) < blast_score_ratio:
                        low_bsr_lines.append(new_line+'\n')
                        apppended += 1
        
        all_matches_lines.append('\n')
        if apppended > 0:
            low_bsr_lines.append('\n')
    
    return [all_matches_lines, low_bsr_lines]


def write_list_to_file(lines, output_file):
    """ Writes a list of strings to a file.
    
        Args:
            lines (list): list with strings as elements, each string representing 
            a line to write to the file.
            output_file (str): complete path to the output file that will be created.
        
        Returns:
            Writes 'output_file' with the strings in the 'lines' list.
    """
    
    with open(output_file, 'w') as out:
        out.writelines(lines)


def allele_set(schema_dictionary):
    """ Retrieves the number of alleles per locus and the identifiers of those alleles.

        Args:
            schema_dictionary (dict): dictionary with locus identifiers as keys and
            with a dictionary with allele identifiers and protein sequences as values.

        Returns:
            alleles_numbers (dict): dictionary with locus identifiers as keys and
            with a list with the total number of alleles and the identifiers of the alleles
            of that locus.
    """

    alleles_numbers = {}
    for locus in schema_dictionary:
        # get alleles identifiers
        alleles = list(schema_dictionary[locus].keys())
        # get total number of alleles
        number_of_alleles = len(alleles)

        alleles_numbers[locus] = [number_of_alleles, alleles]

    return alleles_numbers


def validate_same_locus(matches, alleles_info, blast_score_ratio, bsr_index, genes_ids):
    """ Determines the alleles of the same locus that matched with a representative and the ones that
        did not.

        Args:
            matches (list): list with a dictionary for each locus and that contains the information
            about the BLAST hit between each representative allele of the locus and the hits found
            in the database.
            alleles_info (dict): dictionary with locus identifiers as keys and a list with the
            total number of alleles and the identifiers of all alleles for that locus as values.
            blast_score_ratio (float): the BLAST Score Ratio value used as threshold to
            consider that both sequences in a match belong to the same locus.
            bsr_index (int): index of the BLAST Score Ratio value in the information list
            for every match.

        Returns:
            matched_representatives (dict): dictionary with loci identifiers as keys and
            lists with the information about the number of representative alleles, total number
            of alleles, the number of alleles that matched with a representative allele, if all
            alleles had a match with a representative allele and the identifiers of the alleles that
            had no match with a representative allele.
    """

    # determine passed
    matched_representative = {}
    not_matched = 0
    for dictionary in matches:
        for locus in dictionary:
            # get representative alleles and number of representatives
            representatives = dictionary[locus]
            number_representatives = len(representatives)
            alleles_set = []
            for rep in representatives:
                # get hits for the current representative allele
                hits = representatives[rep]

                for h in hits:
                    # add allele identifier to list if the match has BSR >= 0.6
                    if float(hits[h][bsr_index]) >= blast_score_ratio:
                        alleles_set.append(h)
            
            # determine unique set of allele identifiers with matches with BSR >= 0.6
            alleles_set = set(alleles_set)
            # if the number of allele identifiers that were matched with BSR >= 0.6
            # is equal to the total number of alleles of the locus, it is considered 
            # that all alleles have been validly added to the locus file
            real_allele_num = len(alleles_info[int(locus)][1])
            locus_id = genes_ids[int(locus)]
            if len(alleles_set) == real_allele_num:
                
                matched_representative[locus_id] = [str(number_representatives), str(real_allele_num), 
                                                 str(len(alleles_set))+'/'+str(real_allele_num), 
                                                 'MATCHED', 'NONE']

            # if the number of alleles in both sets differ, some alleles did not have 
            # a match with BSR >= 0.6 with a representative
            elif len(alleles_set) != real_allele_num:

                not_passed = list(set(alleles_info[int(locus)][1]) - set(alleles_set))
                matched_representative[locus_id] = [str(number_representatives), str(real_allele_num), 
                                     str(len(alleles_set))+'/'+str(real_allele_num), 
                                     'NOT MATCHED', not_passed]
                not_matched += 1

    return [matched_representative, not_matched]


def validate_different_loci(matches, representatives, header, blast_score_ratio, genes_ids, alleles_ids):
    """ Creates lines for all matches found between alleles of different loci, for
        all matches between alleles of different loci that had a BSR > blast_score_ratio and
        lines indicating if there were any matches between representative alleles of
        different loci.

        Args:
            matches (dict): dictionary with locus identifiers as keys and a list
            with sublists of matches between alleles of that locus and alleles of another
            locus.
            representatives (list): list with the identifiers of all representative alleles.
            header (str): string that is the columns identifiers for the file.
            blast_score_ratio (float): BLAST Score Ratio value used as threshold to decide
            if matches between alleles of different loci are relevant.

        Returns:
            all_matches_lines (list): list with strings that represent all matches lines to write to a
            file.
            high_bsr_lines (list): list with strings that represent only the matches between alleles
            of different loci that had a BSR > blast_Score_ratio.
            passed (dict): dictionary with loci identifiers as keys and a string
            with the identifiers of matched representatives from other loci (if there
            are any) and if the locus passed the criteria.
    """

    # append header
    all_matches_lines = []
    all_matches_lines.append(header)

    high_bsr_lines = []
    high_bsr_lines.append(header)

    passed = {}
    for d in range(len(matches)):
        # get dictionary for a locus
        all_matches = matches[d]
        # get all matches lists for the current locus
        for locus in all_matches:
            locus_id = genes_ids[int(locus)]
            passed_locus = []
            full_matches = all_matches[locus]
            apppended = 0
            if len(full_matches) > 0:
                for match in full_matches:
                    # create line to add to the list with all matches between
                    # alleles of different loci
                    query_id = alleles_ids[match[0]]
                    hit_id = alleles_ids[match[1]]

                    new_line = ','.join([query_id, hit_id]+match[2:])
                    all_matches_lines.append(new_line+'\n')

                    # append to list if match has BSR >= threshold
                    if float(match[11]) >= blast_score_ratio:
                        high_bsr_lines.append(new_line+'\n')
                        apppended += 1
                        if match[0] and match[1] in representatives:
                            if match[-1] == '11':
                                passed_locus.append(hit_id)

        # if any match between representative alleles of different loci has BSR >= threshold
        # signal that the loci might have some similarities
        if len(passed_locus) > 0:
            passed[locus_id] = '|'.join(passed_locus) + ',NOT_PASSED'
        else:
            passed[locus_id] = 'NONE' + ',PASSED'

        if len(full_matches) > 0:
            all_matches_lines.append('\n')
        if apppended > 0:
            high_bsr_lines.append('\n')

    return [all_matches_lines, high_bsr_lines, passed]


def write_matches_validation(matches, output_file, header):
    """ Writes if a locus passed the criteria that establishes that a locus
        representative should not have a match with a representative of another
        locus with a BSR > 0.6.

        Args:
            matches (dict): dictionary with loci identifiers as keys and information
            about the number of representatives, the number of matched representatives
            and if the locus passed the evaluation.
            output_file (str): path/name of output file to write matches.
            header (str): string with the names of the headers for each column.

        Returns:
            Writes file with three columns:
                - Locus: locus identifier;
                - Number_Representatives: number of representative alleles for the locus;
                - Number_alleles: total number of alleles for the locus;
                - Matched_alleles: number of alleles that had a match with at least one
                representative allele (BSR >= 0.6);
                - Validate: if the locus is valid according to the criteria being evaluated;
                - Unmatched_Alleles: alleles that had not match with a representative allele.
    """

    lines = []
    lines.append(header)
    for locus in matches:
        info = matches[locus]
        if info[-1] != 'NONE':
            info[-1] = '|'.join(info[-1])
            new_line = '{0},{1}\n'.format(locus, ','.join(info))
        else:
            new_line = '{0},{1}\n'.format(locus, ','.join(info))
        lines.append(new_line)

    with open(output_file, 'w') as out:
        out.writelines(lines)


def inter_loci_matches(locus_dictionary):
    """ Determines the matches between alleles of different loci from the full
        set of matches.

        Args:
            locus_dictionary (dict): dictionary with the complete information
            about all loci matches.

        Returns:
            locus_matches (dict): dictionary with genes/loci identifiers as keys
            and lists with complete information about matches between the locus
            that is the key and other loci alleles.
    """

    # for each dictionary/locus
    for locus in locus_dictionary:
        # get the BLAST matches for each allele of that locus
        locus_matches = {locus: []}
        alleles = locus_dictionary[locus]
        # get the hits for each allele
        for allele in alleles:
            hits = alleles[allele]
            # determine if match is between alleles of different loci
            # and keep those matches
            for h in hits:
                query_prefix = hits[h][0].split('_')[0]
                subject_prefix = hits[h][1].split('_')[0]

                if query_prefix != subject_prefix:
                    locus_matches[locus].append(hits[h])

    return locus_matches


def function_helper(input_args):
    """ Runs function by passing set of provided inputs and
        captures exceptions raised during function execution.

        Parameters
        ----------
        input_args : list
            List with function inputs and function object to call
            in the last index.

        Returns
        -------
        results : list
            List with the results returned by the function.
            If an exception is raised it returns a list with
            the name of the function and the exception traceback.
    """

    try:
        results = input_args[-1](*input_args[0:-1])
    except Exception as e:
        func_name = (input_args[-1]).__name__
        traceback_lines = traceback.format_exception(etype=type(e), value=e,
                                                     tb=e.__traceback__)
        traceback_text = ''.join(traceback_lines)
        print('Error on {0}:\n{1}\n'.format(func_name, traceback_text))
        results = [func_name, traceback_text]

    return results


def map_async_parallelizer(inputs, function, cpu, callback='extend',
                           chunksize=1, show_progress=False):
    """ Parallelizes function calls by creating several processes
        and distributing inputs.

        Parameters
        ----------
        inputs : list
            List with inputs to process.
        function
            Function to be parallelized.
        cpu : int
            Number of processes to create (based on the
            number of cores).
        callback : str
            Results can be appended, 'append', to the
            list that stores results or the list of results
            can be extended, 'extend'.
        chunksize : int
            Size of input chunks that will be passed to
            each process. The function will create groups
            of inputs with this number of elements.
        show_progress: bool
            True to show a progress bar with the percentage
            of inputs that have been processed, False
            otherwise.

        Returns
        -------
        results : list
            List with the results returned for each function
            call.
    """

    results = []
    pool = Pool(cpu)
    if callback == 'extend':
        rawr = pool.map_async(function, inputs,
                              callback=results.extend, chunksize=chunksize)
    elif callback == 'append':
        rawr = pool.map_async(function, inputs,
                              callback=results.append, chunksize=chunksize)

    if show_progress is True:
        completed = False
        while completed is False:
            completed = progress_bar(rawr, len(inputs))

    rawr.wait()

    return results


def progress_bar(process, total, tickval=5, ticknum=20, completed=False):
    """ Creates and prints progress bar to stdout.

        Parameters
        ----------
        process : multiprocessing.pool.MapResult
            Multiprocessing object.
        total : int
            Total number of inputs that have to be processed.
        tickval : int
            Progress completion percentage value for each
            tick.
        ticknum : int
            Total number of ticks in progress bar.
        completed : bool
            Boolean indicating if process has completed.

        Returns
        -------
        completed : bool
            Boolean indicating if process has completed.
    """

    # check if process has finished
    if (process.ready()):
        # print full progress bar and satisfy stopping condition
        progress_bar = '[{0}] 100%'.format('='*ticknum)
        completed = True

    # check how many inputs have been processed
    remaining = process._number_left
    if remaining == total:
        # print empty progress bar
        progress_bar = '[{0}] 0%'.format(' '*ticknum)
    else:
        # print progress bar, incremented by 5%
        progress = int(100-(remaining/total)*100)
        progress_tick = progress//tickval
        progress_bar = '[{0}{1}] {2}%'.format('='*progress_tick,
                                              ' '*(ticknum-progress_tick),
                                              progress)

    print('\r', progress_bar, end='')
    time.sleep(0.5)

    return completed


def concatenate_files(files, output_file, header=None):
    """ Concatenates the contents of a set of files.

        Parameters
        ----------
        files : list
            List with the paths to the files to concatenate.
        output_file : str
            Path to the output file that will store the
            concatenation of input files.
        header : str or NoneType
            Specify a header that should be written as the
            first line in the output file.

        Returns
        -------
        output_file : str
            Path to the output file that was created with
            the concatenation of input files.
    """

    with open(output_file, 'w') as of:
        if header is not None:
            of.write(header)
        for f in files:
            with open(f, 'r') as fd:
                shutil.copyfileobj(fd, of)

    return output_file


def flatten_list(list_to_flatten):
    """ Flattens one level of a nested list.

        Parameters
        ----------
        list_to_flatten : list
            List with nested lists.

        Returns
        -------
        flattened_list : str
            Input list flattened by one level.
    """

    flattened_list = list(itertools.chain(*list_to_flatten))

    return flattened_list
