#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Purpose
-------
This module selects representative alleles for a set of loci.

Code documentation
------------------
"""

import os
import sys
import csv
import time
import shutil
import hashlib
import pathlib
import argparse
import traceback
import itertools
import subprocess
import multiprocessing

from Bio import SeqIO
from Bio.Seq import Seq


def file_basename(file_path, file_extension=True):
	"""Extract file basename from a path.

	Parameters
	----------
	file_path : str
		Path to the file.
	file_extension : bool
		Specify if the basename should include the file
		extension.

	Returns
	-------
	basename : str
		File basename extracted from input path.
	"""
	basename = os.path.basename(file_path)

	if file_extension is False:
		# Remove suffix after last '.'
		basename = str((pathlib.Path(basename)).with_suffix(''))

	return basename


def create_directory(directory_path):
	"""Create a diretory if it does not exist."""
	if not os.path.exists(directory_path):
		os.makedirs(directory_path)
		return True
	else:
		return False


def join_paths(parent_path, child_paths):
	"""Create path by joining a parent directory and a list of child paths."""
	joined_paths = os.path.join(parent_path, *child_paths)

	return joined_paths


def read_lines(input_file, strip=True, num_lines=None):
	"""Read lines in a file.

	Parameters
	----------
	input_file : str
		Path to the input file.
	strip : bool
		Specify if lines should be stripped of leading
		and trailing white spaces and new line characters.

	Returns
	-------
	lines : list
		List with the lines read from the input file.
	"""
	with open(input_file, 'r') as infile:
		if num_lines is None:
			lines = [line for line in infile.readlines()]
		else:
			lines = list(islice(infile, num_lines))

	if strip is True:
		lines = [line.strip() for line in lines]

	return lines


def sequence_generator(input_file):
	"""Create a SeqRecord iterator.

	Parameters
	----------
	input_file : str
		Path to a Fasta file.

	Returns
	-------
	records : Bio.SeqIO.FastaIO.FastaIterator
		SeqRecord iterator.
	"""
	# Useful to create the generator
	# Need to exhaust the generator to avoid high memory usage
	records = SeqIO.parse(input_file, 'fasta')

	return records


def hash_sequence(input_string, hash_type='sha256'):
	"""Compute hash of an input string.

	Parameters
	----------
	input_string : str
		Input string to hash.
	hash_type : str
		Hash type/function that will be used to compute the
		hash (any of the hash functions available in the
		hashlib module).

	Returns
	-------
	hashed_string : str
		String representation of the HASH object
		in hexadecimal digits.
	"""
	# get hash function object from hashlib
	hashing_function = getattr(hashlib, hash_type)

	# default encoding is UTF-8
	hashed_string = hashing_function(input_string.encode()).hexdigest()

	return hashed_string


def sequence_lengths(fasta_file, hashed=False):
	"""Determine length of sequences in a FASTA file.

	Read Fasta file and create dictionary with mapping
	between sequence identifiers and sequence lengths.

	Parameters
	----------
	fasta_file : str
		Path to a FASTA file.
	hashed : bool
		If False, sequence headers are used as
		keys. If True, sequence hashes will be
		used as keys.

	Returns
	-------
	lengths : dict
		Dictionary with sequence identifiers as keys and
		sequence lengths as values.
	"""
	records = sequence_generator(fasta_file)
	if hashed is False:
		lengths = {rec.id: len(rec.seq) for rec in records}
	else:
		lengths = {hash_sequence(str(rec.seq)): len(rec.seq) for rec in records}

	return lengths


def fasta_stats(fasta_file):
	"""Determine the number of sequences in a FASTA file and length stats.

	Parameters
	----------
	fasta_file : str
		Path to a FASTA file.

	Returns
	-------
	fasta_file : str
		Path to the FASTA file.
	total_seqs: int
		Total number of records in the FASTA file.
	mean_length: float
		Mean sequence length.
	"""
	seq_lengths = sequence_lengths(fasta_file)
	min_length = min(seq_lengths.values())
	max_length = max(seq_lengths.values())
	mean_length = sum(seq_lengths.values())/len(seq_lengths)
	total_seqs = len(seq_lengths)

	return [fasta_file, total_seqs, min_length, max_length, mean_length]


def distribute_loci(inputs, cores, method):
	"""Create balanced lists of loci to efficiently parallelize function calls.

	Creates balanced lists of loci to distribute per number of
	available cores. Loci lists can be created based on the number
	of sequences per locus (seqcount), the mean length of the
	sequences (length) in each locus or the product of both values
	(seqcount+length).

	Parameters
	----------
	inputs : list
		List with one sublist per locus. Each sublist has
		a locus identifier, the total number of sequences
		and sequence mean legth for that locus.
	cores : int
		The number of loci groups that should be created.
		Based on the number of CPU cores that will be
		used to process the inputs.
	method : str
		"seqcount" to create loci lists based on the total
		number of sequences, "length" to split based
		on mean length of sequences and "seqcount+length" to
		split based on both criteria.

	Returns
	-------
	splitted_ids : list
		List with sublists that contain loci identifiers.
		Sublists are balanced based on the chosen method.
	"""
	# initialize list with sublists to store inputs
	splitted_ids = [[] for cpu in range(cores)]
	# initialize list with chosen criterion values
	# for each sublist of inputs
	splitted_values = [0 for cpu in range(cores)]
	i = 0
	for locus in inputs:
		if method == 'seqcount':
			splitted_values[i] += locus[1]
		elif method == 'length':
			splitted_values[i] += locus[4]
		elif method == 'seqcount+length':
			splitted_values[i] += locus[1] * locus[4]
		splitted_ids[i].append(locus[0])
		# at the end of each iteration, choose the sublist
		# with lowest criterion value
		i = splitted_values.index(min(splitted_values))

	return splitted_ids


def function_helper(input_args):
	"""Run function with provided inputs and capture exceptions.

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
		traceback_lines = traceback.format_exc()
		traceback_text = ''.join(traceback_lines)
		print('\nError on {0}:\n{1}\n'.format(func_name, traceback_text), flush=True)
		results = [func_name, traceback_text]

	return results


def progress_bar(remaining, total, previous, tickval=5, ticknum=20):
	"""Create and print a progress bar to the stdout.

	Parameters
	----------
	remaining : int
		Number of remaining tasks to complete.
	total : int
		Total number of inputs that have to be processed.
	previous : int
		Percentage of tasks that had been completed in the
		previous function call.
	tickval : int
		Progress completion percentage value for each
		tick.
	ticknum : int
		Total number of ticks in progress bar.

	Returns
	-------
	completed : bool
		Boolean indicating if all inputs have been processed.
	"""
	# determine percentage of processed inputs
	progress = int(100-(remaining/total)*100)
	# only print if percentage has changed
	if progress != previous:
		progress_tick = progress//tickval
		progress_bar = '[{0}{1}] {2}%'.format('='*progress_tick,
											  ' '*(ticknum-progress_tick),
											  progress)
		print('\r', progress_bar, end='')

	time.sleep(0.1)

	return progress


def map_async_parallelizer(inputs, function, cpu, callback='extend',
						   chunksize=1, show_progress=False, pool_type='pool'):
	"""Run function in parallel.

	Parameters
	----------
	inputs : list
		List with inputs to process.
	function : func
		Function to be parallelized.
	cpu : int
		Number of processes to create (based on the
		number of CPU cores).
	callback : str
		Results can be appended, "append", to the
		list that stores results or the list of results
		can be extended, "extend".
	chunksize : int
		Size of input chunks that will be passed to
		each process. The function will create groups
		of inputs with this number of elements.
	show_progress : bool
		True to show a progress bar with the percentage
		of inputs that have been processed, False
		otherwise.
	pool_type : str
		The multiprocessing.pool object that will be used,
		Pool or ThreadPool.

	Returns
	-------
	results : list
		List with the results returned for each function
		call.
	"""
	if pool_type == 'pool':
		multiprocessing_function = multiprocessing.pool.Pool
	# Gene prediction uses ThreadPool because Pyrodigal might hang with Pool
	elif pool_type == 'threadpool':
		multiprocessing_function = multiprocessing.pool.ThreadPool

	results = []
	# Use context manager to join and close pool automatically
	with multiprocessing_function(cpu) as pool:
		if callback == 'extend':
			rawr = pool.map_async(function, inputs,
								  callback=results.extend, chunksize=chunksize)
		elif callback == 'append':
			rawr = pool.map_async(function, inputs,
								  callback=results.append, chunksize=chunksize)

		if show_progress is True:
			progress = None
			while progress != 100:
				progress = progress_bar(rawr._number_left, len(inputs), progress)

		rawr.wait()

	return results


def import_sequences(input_file):
	"""Import sequences from a FASTA file.

	Parameters
	----------
	input_file : str
		Path to a FASTA file.

	Returns
	-------
	records_dict : dict
		Dictionary with sequence identifiers as keys and
		sequences as values.
	"""
	records = sequence_generator(input_file)
	# Only want record identifier and sequence, no need to use SeqIO.to_dict
	records_dict = {rec.id: str(rec.seq.upper()) for rec in records}

	return records_dict


def translate_sequence(dna_str, table_id):
	"""Translate a DNA sequence using the BioPython package.

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


def determine_duplicated_seqs(sequences):
	"""Create mapping between sequences and sequence identifiers.

	Parameters
	----------
	sequences : dict
		Dictionary with sequence identifiers as keys and
		sequences as values.

	Returns
	-------
	equal_seqs : dict
		Dictionary with sequences as keys and sequence
		identifiers that are associated with each
		sequence as values.
	"""
	equal_seqs = {}
	for seqid, seq in sequences.items():
		# if protein sequence was already added as key
		if seq in equal_seqs:
			# append new protid
			equal_seqs[seq].append(seqid)
		# else add new protein sequence as key and protid
		# as value
		else:
			equal_seqs[seq] = [seqid]

	return equal_seqs


def determine_longest(seqids, sequences):
	"""Find the longest sequence in a set of sequences.

	Parameters
	----------
	seqids : list
		List with sequence identifiers.
	sequences : dict
		Dictionary with sequence identifiers as keys
		and sequences as values.

	Returns
	-------
	chosen : str
		Sequence identifier of the longest sequence.
	"""
	seqids_tups = [(seqid, sequences[seqid]) for seqid in seqids]
	sorted_tups = sorted(seqids_tups, key=lambda x: len(x[1]), reverse=True)
	chosen = sorted_tups[0][0]

	return chosen


def write_lines(lines, output_file, joiner='\n', write_mode='w'):
	"""Write a list of strings to a file.

	Parameters
	----------
	lines : list
		List with the lines/strings to write to the output
		file.
	output_file : str
		Path to the output file.
	joiner : str
		Character used to join lines.
	write_mode : str
		Specify write mode ('w' creates file if it does not
		exist and truncates and over-writes existing file,
		'a' creates file if it does not exist and appends to
		the end of file if it exists).
	"""
	joined_lines = join_list(lines, joiner)

	write_to_file(joined_lines, output_file, write_mode, '\n')


def make_blast_db(makeblastdb_path, input_fasta, output_path, db_type):
	"""Create a BLAST database.

	Parameters
	----------
	makeblastdb_path : str
		Path to the 'makeblastdb' executable.
	input_fasta : str
		Path to the FASTA file that contains the sequences that
		will be added to the BLAST database.
	output_path : str
		Path to the directory where the database files will be
		created. Database files will have the same basename as
		the `input_fasta`.
	db_type : str
		Type of the database, nucleotide (nuc) or protein (prot).

	Returns
	-------
	stdout : bytes
		BLAST stdout.
	stderr : bytes or str
		BLAST stderr.
	"""
	# Use '-parse-seqids' to be able to specify sequences to align against
	# Use v5 databases (text file with list of sequence IDs needs to be converted with blastdb_aliastool)
	# Decent performance with all BLAST versions, except v2.11 which runs much slower for unkown reasons
	# BLAST <= 2.11 cannot create v4 databases if sequence IDs are alphanumeric and composed of 4 chars
	# v5 databases accept those IDs but replace '-' with '_', which is an issue when chewie is looking for the original IDs
	makedb_cmd = [makeblastdb_path, '-in', input_fasta,
				  '-out', output_path, '-parse_seqids',
				  '-dbtype', db_type, '-blastdb_version', '5']

	makedb_process = subprocess.Popen(makedb_cmd,
									  stdout=subprocess.PIPE,
									  stderr=subprocess.PIPE)

	stdout, stderr = makedb_process.communicate()

	# Exit if it is not possible to create BLAST db
	if len(stderr) > 0:
		sys.exit(f'Could not create BLAST database for {input_fasta}\n'
				 f'{makeblastdb_path} returned the following stderr:\n{stderr}')

	return [stdout, stderr]


def determine_blast_task(sequences, blast_type='blastp'):
	"""Determine the type of BLAST task to execute.

	It is necessary to define the BLAST task if any of the
	sequences to align is shorter that 50 base pairs for
	BLASTn or 30 amino acids for BLASTp.

	Parameters
	----------
	sequences : list
		List that contains strings representing DNA or
		protein sequences.
	blast_type : str
		Used to define the type of application, 'blastn'
		or 'blastp'.

	Returns
	-------
	blast_task : str
		A string that indicates the type of BLAST task to
		execute based on the minimum sequence size.

	Notes
	-----
	More information about the task option at:
		https://www.ncbi.nlm.nih.gov/books/NBK569839/
	"""
	# Get sequence length threshold for BLAST application
	length_threshold = 50 if blast_type == 'blastn' else 30
	sequence_lengths = [len(p) for p in sequences]
	minimum_length = min(sequence_lengths)
	if minimum_length < length_threshold:
		blast_task = '{0}-short'.format(blast_type)
	else:
		blast_task = blast_type

	return blast_task


def run_blast(blast_path, blast_db, fasta_file, blast_output,
			  max_hsps=1, threads=1, ids_file=None, blast_task=None,
			  max_targets=None, composition_stats=None):
	"""Execute BLAST to align sequences against a BLAST database.

	Parameters
	----------
	blast_path : str
		Path to the BLAST application executable.
	blast_db : str
		Path to the BLAST database.
	fasta_file : str
		Path to the FASTA file with sequences to align against
		the BLAST database.
	blast_output : str
		Path to the file that will be created to store the
		results.
	max_hsps : int
		Maximum number of High Scoring Pairs per pair of aligned
		sequences.
	threads : int
		Number of threads/cores used to run BLAST.
	ids_file : str
		Path to a file with sequence identifiers, one per line.
		Sequences will only be aligned to the sequences in the
		BLAST database that match any of the identifiers in this
		file.
	blast_task : str
		Type of BLAST task.
	max_targets : int
		Maximum number of target/subject sequences to align
		against.
	composition_stats : int
		Specify the composition-based statistics method used
		by BLAST.

	Returns
	-------
	stdout : bytes
		BLAST stdout.
	stderr : bytes or str
		BLAST stderr.
	"""
	# Do not retrieve hits with high probability of occuring by chance
	blast_args = [blast_path, '-db', blast_db, '-query', fasta_file,
				  '-out', blast_output, '-outfmt', '6 qseqid qstart qend qlen sseqid slen score',
				  '-max_hsps', str(max_hsps), '-num_threads', str(threads),
				  '-evalue', '0.001']

	# Add file with list of sequence identifiers to align against
	if ids_file is not None:
		blast_args.extend(['-seqidlist', ids_file])
	# Add type of BLASTp or BLASTn task
	if blast_task is not None:
		blast_args.extend(['-task', blast_task])
	# Add maximum number of target sequences to align against
	if max_targets is not None:
		blast_args.extend(['-max_target_seqs', str(max_targets)])
	if composition_stats is not None:
		blast_args.extend(['-comp_based_stats', str(composition_stats)])

	blast_process = subprocess.Popen(blast_args,
								  stdout=subprocess.PIPE,
								  stderr=subprocess.PIPE)

	stdout, stderr = blast_process.communicate()

	# Exit if it is not possible to create BLAST db
	if len(stderr) > 0:
		sys.exit(f'Error while running BLASTp for {fasta_file}\n'
				 f'{blast_path} returned the following error:\n{stderr}')

	return [stdout, stderr]


def read_tabular(input_file, delimiter='\t'):
	"""Read a tabular (TSV) file.

	Parameters
	----------
	input_file : str
		Path to a tabular file.
	delimiter : str
		Delimiter used to separate file fields.

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


def flatten_list(list_to_flatten):
	"""Flatten one level of a nested list.

	Parameters
	----------
	list_to_flatten : list
		Nested list to flatten.

	Returns
	-------
	flattened_list : str
		Input list flattened by one level.
	"""
	flattened_list = list(itertools.chain(*list_to_flatten))

	return flattened_list


def fasta_str_record(record_template, record_data):
	"""Create the string representation of a FASTA record.

	Parameters
	----------
	record_template : str
		String template to construct the FASTA record.
	record_data : list
		List with the elements to add to the string.

	Returns
	-------
	record : str
		String representation of the FASTA record.
	"""
	record = record_template.format(*record_data)

	return record


def fasta_lines(template, records_data):
	"""Create a list with FASTA records.

	Parameters
	----------
	template : str
		String template to construct the FASTA record.
	records_data : list
		A list with one sublist per FASTA record.
		Each sublist contains the elements to insert
		inside the template placeholders.

	Returns
	-------
	seqs_lines : list
		A list with strings representing FASTA records.
	"""
	seqs_lines = [fasta_str_record(template, arg) for arg in records_data]

	return seqs_lines


def run_blastdb_aliastool(blastdb_aliastool_path, seqid_infile, seqid_outfile):
	"""Convert list of sequence identifiers into binary format.

	Parameters
	----------
	blastdb_aliastool_path : str
		Path to the blastdb_aliastool executable.
	seqid_infile :  str
		Path to the file that contains the list of sequence identifiers.
	seqid_outfile : str
		Path to the output file in binary format to pass to the -seqidlist
		parameter of BLAST>=2.10.

	Returns
	-------
	stdout : bytes
		BLAST stdout.
	stderr : bytes or str
		BLAST stderr.
	"""
	blastdb_aliastool_args = [blastdb_aliastool_path, '-seqid_file_in',
							  seqid_infile, '-seqid_file_out', seqid_outfile]

	blastdb_aliastool_process = subprocess.Popen(blastdb_aliastool_args,
												 stdout=subprocess.PIPE,
												 stderr=subprocess.PIPE)

	stdout, stderr = blastdb_aliastool_process.communicate()

	# Exit if it is not possible to create BLAST db
	if len(stderr) > 0:
		sys.exit(f'Could not convert {seqid_infile} to binary format.\n'
				 f'{blastdb_aliastool_path} returned the following error:\n{stderr}')

	return [stdout, stderr]


def join_list(lst, delimiter):
	"""Join all elements in a list into a single string.

	Parameters
	----------
	lst : list
		List with elements to be joined.
	delimiter : str
		Character used to join list elements.

	Returns
	-------
	joined_list : str
		A single string with all elements in the input
		list joined by the character chosen as link.
	"""
	joined_list = delimiter.join(lst)

	return joined_list


def write_to_file(text, output_file, write_mode, end_char):
	"""Write a single string to a file.

	Parameters
	----------
	text : str
		A single string to write to the output file.
	output_file : str
		Path to the output file.
	write_mode : str
		Specify write mode ('w' creates file if it does not
		exist and truncates and over-writes existing file,
		'a' creates file if it does not exist and appends to
		the end of file if it exists.).
	end_char : str
		Character added to the end of the file.
	"""
	with open(output_file, write_mode) as out:
		out.write(text+end_char)


def bsr_categorizer(blast_results, representatives,
					representatives_scores, min_bsr, max_bsr):
	"""Identify BLAST hits below and above the BSR min and max thresholds.

	Parameters
	----------
	blast_results : list of list
		A list with sublists, each sublist contains information
		about a BLAST hit.
	representatives : list
		List with sequence identifiers of representative
		sequences.
	representatives_scores : dict
		Dictionary with self BLAST raw score for every
		representative.
	min_bsr : float
		Minimum BSR value accepted to consider a sequence as
		a possible new representative.
	max_bsr : float
		Maximum BSR value accepted to consider a sequence as
		a possible new representative.

	Returns
	-------
	List with the following elements:
		high_bsr : list
			list with all sequence identifiers of subject
			sequences that had hits with a BSR higher than the
			maximum defined threshold.
		low_bsr : list
			list with all sequence identifiers of subject
			sequences that had hits with a BSR lower than the
			minimum defined threshold.
	"""
	high_bsr = []
	hotspot_bsr = []
	low_bsr = []

	high_reps = {}
	hot_reps = {}
	low_reps = {}

	filtered_results = [res for res in blast_results
						if res[0] != res[4] and res[4] not in representatives]
	bsr_values = [float(res[-1])/representatives_scores[res[0]]
				  for res in filtered_results]

	high_bsr = [res[4] for ind, res in enumerate(filtered_results)
				if bsr_values[ind] >= max_bsr]
	low_bsr = [res[4] for ind, res in enumerate(filtered_results)
			   if bsr_values[ind] < min_bsr]
	hotspot_bsr = [res[4] for ind, res in enumerate(filtered_results)
				   if bsr_values[ind] >= min_bsr and bsr_values[ind] < max_bsr]

	for ind, res in enumerate(filtered_results):
		if bsr_values[ind] >= min_bsr:
			high_reps.setdefault(res[0], []).append(res[4])
		if bsr_values[ind] < min_bsr:
			low_reps.setdefault(res[0], []).append(res[4])
		if bsr_values[ind] >= min_bsr and bsr_values[ind] < max_bsr:
			hot_reps.setdefault(res[0], []).append(res[4])

	# Identify representatives that only led to low BSR
	low_reps = list(set(low_reps) - set(high_reps))

	return [high_bsr, low_bsr, hotspot_bsr, high_reps, low_reps, hot_reps]


def select_candidate(candidates, proteins, seqids,
					 representatives, final_representatives):
	"""Select a new representative sequence.

	Parameters
	----------
	candidates : list
		List with the sequence identifiers of all candidates.
	proteins : dict
		A dictionary with protein identifiers as keys and
		protein sequences as values.
	seqids : list
		A list with the sequence identifiers that still have
		no representative (representatives identifiers are
		included because they have to be BLASTed in order to
		determine their self score).
	representatives : list
		The sequence identifiers of all representatives.

	Returns
	-------
	representatives : list
		The set of all representatives, including the new
		representative that was chosen by the function.
	"""
	# With more than one sequence as candidate, select longest
	if len(candidates) > 1:
		# Determine length of all candidates
		candidates_len = [(seqid, len(proteins[seqid]))
						  for seqid in candidates]

		# Order representative candidates by length descending order
		candidates_len = sorted(candidates_len, key=lambda x: x[1],
								reverse=True)

		# Longest allele is the new representative
		representatives.append(candidates_len[0][0])
		final_representatives.append(candidates_len[0][0])
	# If there is only one candidate, keep that
	elif len(candidates) == 1:
		representatives.append(candidates[0])
		final_representatives.append(candidates[0])
	# If no hit qualifies and there are still sequences
	# without representative
	elif len(candidates) == 0 and \
			len(seqids) > len(representatives):
		# Determine length of remaining sequences
		# (representatives not included)
		candidates_len = [(seqid, len(proteins[seqid]))
						  for seqid in seqids
						  if seqid not in representatives]
		# Sort by descending length
		candidates_len = sorted(candidates_len, key=lambda x: x[1],
								reverse=True)
		# Longest of remaining sequences is new representative
		representatives.append(candidates_len[0][0])
		final_representatives.append(candidates_len[0][0])

	return [representatives, final_representatives]


def adapt_loci(loci, schema_path, schema_short_path, bsr,
			   table_id, blastp_path, makeblastdb_path,
			   blastdb_aliastool_path):
	"""Adapts a set of loci from an external schema.

	Adapts an external schema for usage with chewBBACA. Removes invalid
	alleles and selects representative alleles to include in the "short"
	directory.

	Parameters
	----------
	loci_list : list
		A list with the following elements:

		- List with paths to the files to be processed.
		- Path to the schema directory.
		- Path to the "short" directory.
		- BLAST Score Ratio value.
		- Genetic code.

	Returns
	-------
	The function writes the schema files .
	"""
	for locus in loci:
		representatives = []
		final_representatives = []
		rep_self_scores = {}

		# Get locus identifier (does not include extension)
		locus_id = file_basename(locus, file_extension=False)

		# Create paths to gene files in new schema
		locus_file = join_paths(schema_path, [f'{locus_id}.fasta'])
		locus_short_file = join_paths(schema_short_path, [f'{locus_id}_short.fasta'])

		# Create temp directory for current gene
		locus_temp_dir = join_paths(schema_path, [f'{locus_id}_temp'])
		# Create temp directory for the current gene
		create_directory(locus_temp_dir)

		# Dictionaries mapping gene identifiers to DNA and Protein sequences
		locus_seqs = import_sequences(locus)
		prot_seqs = {k: str(translate_sequence(v, table_id)) for k, v in locus_seqs.items()}

		if len(locus_seqs) > 1:
			# Identify DNA sequences that code for same protein
			equal_prots = determine_duplicated_seqs(prot_seqs)

			# Get only one identifier per protein
			ids_to_blast = [protids[0] for protein, protids in equal_prots.items()]

			# Get longest sequence as first representative
			longest = determine_longest(ids_to_blast, prot_seqs)
			representatives.append(longest)
			final_representatives.append(longest)

			# Create FASTA file with distinct protein sequences
			protein_file = join_paths(locus_temp_dir, [f'{locus_id}_protein.fasta'])
			protein_data = [[i, prot_seqs[i]] for i in ids_to_blast]
			protein_lines = fasta_lines('>{0}\n{1}', protein_data)
			write_lines(protein_lines, protein_file)

			# Create blastdb with all distinct proteins
			blastp_db = join_paths(locus_temp_dir, [locus_id])
			db_std = make_blast_db(makeblastdb_path, protein_file, blastp_db, 'prot')

			# Determine appropriate blastp task (proteins < 30aa need blastp-short)
			blastp_task = determine_blast_task(equal_prots)

			# Cycle BLAST representatives against non-representatives until
			# all non-representatives have a representative
			while len(set(ids_to_blast) - set(representatives)) != 0:
				# create FASTA file with representative sequences
				rep_file = join_paths(locus_temp_dir, [f'{locus_id}_rep_protein.fasta'])
				rep_protein_data = [[r, prot_seqs[r]] for r in representatives]
				rep_protein_lines = fasta_lines('>{0}\n{1}', rep_protein_data)
				write_lines(rep_protein_lines, rep_file)

				# Compute self-score for representative alleles
				for seqid in representatives:
					if seqid not in rep_self_scores:
						record = fasta_str_record('>{0}\n{1}', [seqid, prot_seqs[seqid]])
						current_rep_file = join_paths(locus_temp_dir, [f'{seqid}_solo.fasta'])
						write_lines([record], current_rep_file)
						# Create file with representative seqid to only compare against self
						id_file = join_paths(locus_temp_dir, [f'{seqid}_ids.txt'])
						write_lines([seqid], id_file)
						binary_file = f'{id_file}.bin'
						blast_std = run_blastdb_aliastool(blastdb_aliastool_path, id_file, binary_file)
						id_file = binary_file

						rep_blastout = join_paths(locus_temp_dir, [f'{seqid}_blastout.tsv'])
						# Cannot get self-alignemnt for some sequences if composition-based stats is enabled
						blast_std = run_blast(blastp_path, blastp_db, current_rep_file,
											  rep_blastout, 1, 1,
											  id_file, 'blastp', None, 0)
						rep_results = read_tabular(rep_blastout)
						if len(rep_results) > 0:
							rep_self_scores[rep_results[0][0]] = float(rep_results[0][6])
						else:
							print('Could not determine the self-alignment raw '
								f'score for {rep_results[0][0]}')

				# Create file with seqids to BLAST against
				ids_str = join_list([str(i) for i in ids_to_blast if i not in representatives], '\n')
				ids_file = join_paths(locus_temp_dir, [f'{locus_id}_ids.txt'])
				write_to_file(ids_str, ids_file, 'w', '')
				binary_file = f'{ids_file}.bin'
				blast_std = run_blastdb_aliastool(blastdb_aliastool_path, ids_file, binary_file)
				ids_file = binary_file

				# BLAST representatives against non-represented
				blast_output = join_paths(locus_temp_dir, [f'{locus_id}_blast_out.tsv'])
				# Set 'max_target_seqs' to huge number because BLAST only
				# returns 500 hits by default
				blast_std = run_blast(blastp_path, blastp_db, rep_file,
									  blast_output, 1, 1, ids_file,
									  blastp_task, 100000)

				# Import BLAST results
				blast_results = read_tabular(blast_output)

				# Divide results into high, low and hot BSR values
				hitting_high, hitting_low, hotspots, high_reps, low_reps, hot_reps = \
					bsr_categorizer(blast_results, representatives,
									rep_self_scores, bsr, bsr+0.1)

				excluded_reps = []
				# Remove high BSR hits that have representative
				hitting_high = set(hitting_high)
				ids_to_blast = [i for i in ids_to_blast if i not in hitting_high]

				# Remove representatives that led to high BSR with subjects that were removed
				prunned_high_reps = {k: [r for r in v if r in ids_to_blast] for k, v in high_reps.items()}
				reps_to_remove = [k for k, v in prunned_high_reps.items() if len(v) == 0]

				excluded_reps.extend(reps_to_remove)

				# Determine smallest set of representatives that allow to get all cycle candidates
				excluded = []
				hotspot_reps = set(flatten_list(list(hot_reps.values())))
				for rep, hits in hot_reps.items():
					common = hotspot_reps.intersection(set(hits))
					if len(common) > 0:
						hotspot_reps = hotspot_reps - common
					else:
						excluded.append(rep)

				excluded_reps.extend(excluded)

				# Remove representatives that only led to low BSR
				excluded_reps.extend(low_reps)

				representatives = [rep for rep in representatives if rep not in excluded_reps]
				ids_to_blast = [i for i in ids_to_blast if i not in excluded_reps]

				# Determine next representative from candidates
				rep_candidates = list(set(hotspots) - hitting_high)
				# Sort to guarantee reproducible results with same datasets
				rep_candidates = sorted(rep_candidates, key=lambda x: int(x))
				representatives, final_representatives = select_candidate(rep_candidates,
																		  prot_seqs,
																		  ids_to_blast,
																		  representatives,
																		  final_representatives)

				# Remove files created for current gene iteration
				os.remove(rep_file)
				os.remove(blast_output)
				os.remove(ids_file)
		else:
			final_representatives = list(prot_seqs.keys())

		# Write schema file with all alleles
		locus_data = [[k, v] for k, v in locus_seqs.items()]
		locus_lines = fasta_lines('>{0}\n{1}', locus_data)
		write_lines(locus_lines, locus_file)

		# Get total number of valid sequences
		total_sequences = len(locus_lines)

		# Write schema file with representatives
		locus_rep_data = [[r, locus_seqs[r]]
						  for r in final_representatives]
		locus_rep_lines = fasta_lines('>{0}\n{1}', locus_rep_data)
		write_lines(locus_rep_lines, locus_short_file)

		# Get number of representatives
		representatives_number = len(locus_rep_lines)

		shutil.rmtree(locus_temp_dir)

	return True


def main(input_file, output_directory, cpu_cores, blast_score_ratio,
		 translation_table):
	"""
	Adapt a schema to be used with chewBBACA.

	Parameters
	----------
	input_file : str
		Path to a TXT file with the list of schema loci to adapt.
	output_directory :  list
		Path to the output directory to create (the main schema
		directory and the 'short' directory to store representative
		alleles).
	cpu_cores : int
		Number of CPU cores that will be used to run the process.
	blast_score_ratio : float
		The BLAST Score Ratio value that will be used to evaluate
		allele similarity and select representative alleles.
	translation_table : int
		Genetic code used to translate alleles.
	"""
	# Define output paths
	schema_path = os.path.abspath(output_directory)
	schema_short_path = join_paths(schema_path, ['short'])

	# Create output directories
	create_directory(schema_path)
	create_directory(schema_short_path)

	# Import list of loci to adapt
	loci_list = read_lines(input_file, strip=True)
	print(f'Number of loci to adapt: {len(loci_list)}')
	# Count number of sequences and mean length per locus
	loci_info = []
	loci_pools = multiprocessing.Pool(processes=cpu_cores)
	gp = loci_pools.map_async(fasta_stats, loci_list, callback=loci_info.extend)
	gp.wait()

	# Split files according to number of sequences and sequence mean length
	# Divide into 100 input sets to get 1% progress resolution
	even_loci_groups = distribute_loci(loci_info, 100, 'seqcount')
	# With few inputs, some sublists might be empty
	even_loci_groups = [i for i in even_loci_groups if len(i) > 0]
	# Add common arguments
	blastp_path = 'blastp'
	makeblastdb_path = 'makeblastdb'
	blastdb_aliastool_path = 'blastdb_aliastool'
	even_loci_groups = [[i, schema_path, schema_short_path,
						 blast_score_ratio, translation_table,
						 blastp_path, makeblastdb_path, blastdb_aliastool_path,
						 adapt_loci] for i in even_loci_groups]

	print(f'Adapting...')
	adaptation_data = map_async_parallelizer(even_loci_groups,
											 function_helper,
											 cpu_cores,
											 show_progress=True)

	print('\nDone.')

	return True

if __name__ == '__main__':
	pass
