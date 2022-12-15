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
    # use context manager to join and close pool automatically
    with Pool(cpu) as pool:
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
