#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
AUTHOR

    Rafael Mamede
    github: @rfm-targa

DESCRIPTION

    This module has the collection of functions that are used to validate schemas.
    The schema structure should respect a set of criteria that are used during the
    CreateSchema and AlleleCall processes. All alleles for a gene present in the schema
    should be selected only if they share a BSR >= 0.6 with at least one representative
    allele of that gene. Representative alleles from different loci should not share
    a BSR >= 0.6. However, the process does not guarantee that all alleles for a gene 
    in the schema share a BSR >= 0.6 or that non-representative sequences from one gene
    cannot have a BSR >= 0.6 with alleles from other gene.
    The functions in this module were created with the intent of validating schemas and
    help determine if a given schema or subset of a schema respects the necessary criteria
    and does not deviate too much from those criteria due to aspects that are not controlled
    during the CreateSchema or AlleleCall processes.

    Several functions in this module require that 'Biopython' be installed in the working
    environment.
"""


import os
import sys
import csv
import shutil
from collections import Counter

from Bio import SeqIO
from Bio.Seq import Seq


def schema_dictionary(schema_files, files_directory):
    """ Creates a dictionary that stores all allele sequences for each gene of
        the schema.
    
        Args: 
            schema_files (list): the list of files names that constitute the schema.
            files_directory (str): the directory with the schema files.
        
        Returns: 
            schema_dict (dict): a dictionary with gene identifiers as keys. Each 
            key has a sub dictionary with alleles identifiers as keys and alleles 
            sequences as values.
        
        Example:query_prefix
            
            >>> schema_dictionary(['gene_1.fasta','gene_2.fasta'], '/home/user/schema_folder')
            {'gene_1-protein1':{'gene_1-protein1_allele1':'ATGAAATAA',
                                'gene_1-protein1_allele2':'ATGAGGTAA'},
             'gene_1-protein2':{'gene_1-protein2_allele1':'ATGACATAA',
                                'gene_1-protein2_allele2':'ATGAGCTAA'},...}
    """

    schema_dict = {}
    # for each FASTA file that has the alleles for a locus
    gene_ids_dict = {}
    gene_num = 1
    alleles_ids_dict = {}
    for file in schema_files:
        gene_name = file.split('.')[0]
        gene_ids_dict[gene_num] = gene_name
        # create a sub dictionary to hold the allele sequences for each gene
        schema_dict[gene_num] = {}
        file_name = '{0}/{1}'.format(files_directory, file)
        allele_num = 1
        # get the allele identifier and DNA sequence for each allele in the file
        for allele in SeqIO.parse(file_name, 'fasta'):
            allele_id = allele.id
            allele_numid = '{0}_{1}'.format(gene_num, allele_num)
            alleles_ids_dict[allele_numid] = allele_id
            allele_sequence = str(allele.seq)

            schema_dict[gene_num][allele_numid] = allele_sequence
            allele_num += 1

        gene_num += 1

    return [schema_dict, gene_ids_dict, alleles_ids_dict]


def short_schema_dict(schema_files, files_directory, genes_ids, alleles_ids):
    """
    """
    
    # invert genes and alleles identifiers dictionaries
    genes_ids = {v:k for k, v in genes_ids.items()}
    alleles_ids = {v:k for k, v in alleles_ids.items()}
    
    schema_dict = {}
    # for each FASTA file that has the alleles for a locus
    for file in schema_files:
        # need to find a way to deal with all kinds of identifiers...
        gene_name = file.split('.')[0]
        gene_name = gene_name.rstrip('_short')
        gene_num = genes_ids[gene_name]
        # create a sub dictionary to hold the allele sequences for each gene
        schema_dict[gene_num] = {}
        file_name = '{0}/{1}'.format(files_directory, file)
        # get the allele identifier and DNA sequence for each allele in the file
        for allele in SeqIO.parse(file_name, 'fasta'):
            allele_id = alleles_ids[allele.id]
            allele_sequence = str(allele.seq)

            schema_dict[gene_num][allele_id] = allele_sequence
    
    return schema_dict


def translate_schema(schema_dictionary):
    """ Creates a dictionary that stores the protein sequences coded in the DNA 
        sequences of each allele per gene in the schema.

        Args: 
            schema_dictionary (dict): a dictionary with gene identifiers as keys. Each 
            key has a sub dictionary with alleles identifiers as keys and alleles 
            sequences as values.

        Returns: 
            protein_schema_dict (dict): a dictionary with gene identifiers as keys. Each 
            key has a sub dictionary with alleles identifiers as keys and alleles 
            protein sequences as values.

        Example:

            >>> translate_schema({'gene_1-protein1':{'gene_1-protein1_allele1':'ATGAAATAA',
                                                     'gene_1-protein1_allele2':'ATGAGGTAA'},...})
            {'gene_1-protein1':{'gene_1-protein1_allele1':'MK',
                                'gene_1-protein1_allele2':'MR'},...}
    """

    protein_schema_dict = {}
    for gene_name in schema_dictionary:
        # get the dictionary with the DNA sequences of the current gene
        allele_dict = schema_dictionary[gene_name]
        for allele, sequence in allele_dict.items():
            # translate the DNA sequences of all alleles
            translation = translate(sequence)
            # only proceed if the sequence could be translated
            if isinstance(translation, list):
                translated_sequence = translation[0]
                # add protein sequence of allele
                if gene_name in protein_schema_dict:
                    protein_schema_dict[gene_name][allele] = translated_sequence
                # add new gene key to dict and first allele protein sequence
                else:
                    protein_schema_dict[gene_name] = {allele: translated_sequence}
        
    return protein_schema_dict


def schema_fasta_lines(protein_schema_dictionary):
    """ Creates lines with headers and protein sequences typical of a FASTA file.
    
        Args: 
            protein_schema_dictionary (dict): a dictionary with gene identifiers as keys. Each
            key has a sub dictionary with alleles identifiers as keys and alleles protein 
            sequences as values.
        
        Returns:
            protein_lines (dict): a dictionary with locus identifiers as keys and a list 
            of strings with headers and protein sequences as values.
        
        Example:
            
            >>> schema_fasta_lines({'gene_1-protein1':{'gene_1-protein1_allele1':'MK',
                                                     'gene_1-protein1_allele2':'MR'},...})
            {gene_1-protein1: ['> gene_1-protein1_allele1', 'MK',
                               '> gene_1-protein1_allele2', 'MR',...],
            {gene_2-protein2: ['> gene_2-protein2_allele1', 'ML',
                               '> gene_2-protein2_allele2', 'MN',...],
            ...}
    """

    protein_lines = {}
    # for each sequence identifier
    for gene_name in protein_schema_dictionary:
        # get the dictionary with the alleles for that gene
        allele_dict = protein_schema_dictionary[gene_name]
        alleles_lines = []
        # for each allele sequence identifier and allele sequence
        for sequence_id, sequence in allele_dict.items():
            # create the header and sequence lines and append them sequentially
            header = '> {}'.format(sequence_id)
            protein = sequence

            alleles_lines.append(header)
            alleles_lines.append(protein)
        
        protein_lines[gene_name] = alleles_lines

    return protein_lines


def reverse_complement(dna_sequence):
    """ Determines the reverse complement of given DNA strand.
    
        Args:
            dna_sequence (str): string representing a DNA sequence.
        
        Returns:
            reverse_complement_strand (str): the reverse complement of the DNA sequence, without
            lowercase letters.
        
        Example:
            >>> reverse_complement('ATCGgcaNn')
            'NNTGCCGAT'
    """
    
    base_complement = {'A':'T', 'C':'G', 'G':'C', 'T':'A',
                      'a':'T', 'c':'G', 'g':'C', 't':'A',
                      'n':'N', 'N':'N'}
    
    # convert string into list with each character as a separate element
    sequence_bases = list(dna_sequence)
    
    # determine complement strand
    sequence_complement = [base_complement[base] for base in sequence_bases]
    
    complement_strand = ''.join(sequence_complement)
    
    # reverse strand
    reverse_complement_strand = complement_strand[::-1]
    
    return reverse_complement_strand


def translate(dna_sequence):
    """ Converts a given DNA sequence into a protein sequence. The strand will 
        be translated according to the sense strand if it has start and stop codons, 
        otherwise it will be reverse complemented and translated if possible.
        
        Args:
            dna_sequence (str): a string representing a DNA sequence.
        
        Returns:
            protseq (str): a string representing the protein sequence obtained from 
            the translation of the input DNA sequence.
            translated_strand (str): indicates the DNA strand that coded for the protein
            sequence (sense or antisense).
            
        Example:
            >>> translate("GTGACACCAAAACCATGA")
            ['MTPKP', 'sense']
        
        Note:
            The sequence is scanned for ambiguous nucleotides according to IUPAC definitions.
            The translation table used is the 'The Bacterial, Archaeal and Plant Plastid Code'
            from NCBI.
            For more information about the translation table used visit:
                https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi
    """
    
    sequence = dna_sequence
    tableid = 11
    # ATC and ATA are start codons? OMG must change this function asap
    start_codons = ('ATG','GTG','TTG','ATT','CTG','ATC','ATA')
    stop_codons = ('TAA','TAG','TGA')
    translated_strand = ''

    # sequences with ambiguous nucleotides are not translated
    if any((nuc in 'RYWSMKHBVDN') for nuc in sequence):
        ambiguous = 'amb'
        
        return ambiguous
    
    # if the sequence has no ambiguous nucleotides
    else:
        # check for start and stop codons
        start = sequence.startswith(start_codons)
        stop = sequence.endswith(stop_codons)
        
        # translate sequence if it is a CDS
        if start and stop:
            myseq = Seq(sequence)
            protseq = Seq.translate(myseq, table=tableid, cds=True)
            translated_strand = 'sense'
        
        # reverse complement the strand to check for CDS if sense strand had no CDS
        else:
            revC_seq = reverse_complement(sequence)
            start = revC_seq.startswith(start_codons)
            stop = revC_seq.endswith(stop_codons)

            if start and stop:
                myseq = Seq(revC_seq)
                protseq = Seq.translate(myseq, table=tableid, cds=True)
                translated_strand = 'antisense'

    if translated_strand == '':
        return 1
        
    return [str(protseq), translated_strand]


def create_blastdb(write_directory, input_file, verbose=False):
    """ Creates a BLASTp database based on the input sequences file.
        
        Args:
            write_directory (str): path to the directory where the new BLASTp database 
            directory should be created.
            input_file (str): path to the protein sequences file in FASTA format. The 
            basename of this string is used as prefix for the name of the database.
        
        Returns:
            database_name (str): name of the BLASTp database that is created.
            
        Example:
            >>> create_blastdb('/home/user/workdir/', '/home/user/workdir/file.fasta')
            file_db
    """
    
    basename = os.path.splitext(os.path.basename(input_file))[0]
    blastdb_dirname = '{0}/{1}_blastdbs'.format(write_directory, basename)
    database_name = '{0}/{1}_db'.format(blastdb_dirname, basename)
    
    database_files = [database_name + '.pin', database_name + '.phr', database_name + '.psq']

    if not os.path.isdir(blastdb_dirname):
        os.makedirs(blastdb_dirname)
    
    protein_db_cmd = 'makeblastdb -in {0} -out {1} -dbtype prot -logfile {1}_blast.log'.format(input_file, database_name)

    if not all([os.path.isfile(file) for file in database_files]):
        
        os.system(protein_db_cmd)
        
        if not all([os.path.isfile(file) for file in database_files]):
            raise Exception('BLASTp database could not be created, one or more database files could not be created.')
        elif not all([os.stat(file).st_size > 0 for file in database_files]):
            raise Exception('BLASTp database could not be created, one or more database files are empty.')
        elif verbose == True:
            print('BLASTp database succesfully created!')

    else:
        print('BLASTp database files found. Using existing database files..')
        
    return database_name


def write_fasta(fasta_lines, output_file):
    """ Writes strings/lines stored in an input list to a file, sequentially.
        Lines should represent FASTA headers and sequences.
        
        Args:
            fasta_lines (list): a list with lines/strings.
            output_file (str): the name of the output file.
                
        Returns:
            Creates 'output_file' FASTA file and prints the number of sequences
            written into the FASTA file.
    """
    
    joined_lines = '\n'.join(fasta_lines)
    #number_sequences = len(fasta_lines) // 2
    
    with open(output_file, 'w') as file:
        file.write(joined_lines)
    
    # Inform user of how many sequences were written into the file.
    #print('Wrote {0} FASTA sequences to {1}'.format(number_sequences, output_file))


def write_per_locus_fasta_file(loci_lines, output_directory, suffix):
    """ Writes a set of FASTA files based on a dictionary with gene identifiers as 
        keys and a list of headers and sequences for each gene.
        
        Args:
            lines (dict): a dictionary with genes/loci identifiers as keys and a list 
            of headers and sequences as values for each key.
            output_directory (str): the path to the directory where FASTA files should be created.
                
        Returns:
            fasta_paths (list): a list with the full path for each FASTA file created.
    """

    fasta_paths = []
    for locus in loci_lines:
        fasta_lines = loci_lines[locus]
        output_file = '{0}/{1}{2}.fasta'.format(output_directory, locus, suffix)
        fasta_paths.append(output_file)
        write_fasta(fasta_lines, output_file)
    
    # sort paths in list based on the locus name so that the lists of paths
    # for the files with all alleles and for the files with the representative alleles
    # have the same order
    if 'short' in fasta_paths[0]:
        sorting_function = lambda x: (x.split('/')[-1]).split('_short')[0]
    else:
        sorting_function = lambda x: (x.split('/')[-1]).split('.fasta')[0]
    
    fasta_paths.sort(key = sorting_function)
    
    return fasta_paths


def write_all_loci_fasta_file(loci_lines, output_directory):
    """ Joins the headers and sequences lines of all genes/loci and writes 
        them into file.
        
        Args:
            lines (dict): a dictionary with genes/loci identifiers as keys and a list 
            of headers and sequences as values for each key.
            output_directory (str): the path to the directory where FASTA files should be created.
        
        Returns:
            Creates a FASTA file with headers and sequences of all genes/loci.
    """
    
    fasta_lines = []
    for locus in loci_lines:
        fasta_lines += loci_lines[locus]
    
    write_fasta(fasta_lines, output_directory)


def within_loci_blast(blast_input):
    """ Creates a BLASTp database with the full set of alleles in the schema and 
        BLASTs the representative sequences against the BLASTp database.
        
        Args:
            full_fastas (list): a list with the full paths for the FASTA files with all 
            alleles.
            short_fastas (list): a list with the full paths for the FASTA files with the 
            representative sequences.
            output_directory (str): name of the directory where the BLASTp daabase should be 
            created.
            num_threads (int): number o CPU cores to use to run BLASTp.
        Returns:
            blast_files (list): a list with the paths to the BLAST output files created.
    """

    full_file = blast_input[0]
    short_file = blast_input[1]
    output_directory = blast_input[2]
    num_threads = blast_input[3]
    
    # determine if minimum allele length is smaller than 30aa
    small = False
    for allele in SeqIO.parse(full_file, 'fasta'):
        if len(str(allele.seq)) < 30:
            small = True
    
    # Create BLASTp database with all sequences
    blastdb_name = create_blastdb(output_directory, full_file)
    
    filename = os.path.basename(full_file).split('.')[0]
    
    blast_out_file = '{0}/{1}_blast_out.tsv'.format(output_directory, filename)
    
    # do not forget to define max HSPs (High-scoring Segment Pairs) or BLAST might return several hits per subject with the same query
    # increase number of maximum targets because BLAST default is 500
    if small is False:
        blast_command = ('blastp -db {0} -query {1} -out {2} '
                        '-outfmt "6 qseqid sseqid pident length qstart qend sstart send evalue bitscore score" '
                        '-max_hsps 1 -num_threads {3} -max_target_seqs 100000'.format(blastdb_name, short_file, blast_out_file, num_threads))
    # if there are proteins smaller than 30aa, use blastp-short
    elif small is True:
        blast_command = ('blastp -task blastp-short -db {0} -query {1} -out {2} '
                        '-outfmt "6 qseqid sseqid pident length qstart qend sstart send evalue bitscore score" '
                        '-max_hsps 1 -num_threads {3} -max_target_seqs 100000'.format(blastdb_name, short_file, blast_out_file, num_threads))
    
    os.system(blast_command)
    
    os.remove(full_file)
    os.remove(short_file)
    
    blastdb_dirname = os.path.dirname(blastdb_name)
    shutil.rmtree(blastdb_dirname)
    
    return blast_out_file


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


def determine_representatives(representatives_schema):
    """ Retrieves the identifiers of all representative alleles in the input schema.

        Args:
            representative_schema (dict): dictionary with locus identifiers as keys and
            a dictionary with alleles identifiers as keys and sequences as values.

        Returns:
            representative_identifiers (list): list with all representative alleles
            identifiers of all loci.
    """

    representatives_identifiers = []
    # for each allele identifier in the short schema
    for locus in representatives_schema:
        # add each representative allele identifier to the same list
        representatives = list(representatives_schema[locus].keys())
        representatives = [rep.split('.')[0] for rep in representatives]
        representatives_identifiers.extend(representatives)

    return representatives_identifiers


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


def identify_input_type(input_argument):
    """ Determines if the input argument path is for a file or for a directory.

        Args:
            input_argument (str): path to a file or directory.

        Returns:
            A list with two variables:
                input_files (list): list of files names to all schema files to use.
                schema_directory (str): path of the directory with the schema files.
    """

    # if it is a directory, list all files in that directory
    if os.path.isdir(input_argument) is True:
        input_files = [f for f in os.listdir(input_argument) if f.endswith('.fasta') is True]
        schema_directory = input_argument

    # if it is a file, get the basename for all files
    elif os.path.isfile(input_argument) is True:
        with open(input_argument, 'r') as file:
            input_files = [f.strip() for f in file.readlines()]
            schema_directory = os.path.dirname(input_files[0])
            input_files = [os.path.basename(f) for f in input_files]

    return [input_files, schema_directory]
