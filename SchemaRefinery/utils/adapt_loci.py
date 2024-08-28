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
import time
import shutil
import traceback
import multiprocessing
    
from Bio import SeqIO
from Bio.Seq import Seq

try:
    from utils import (file_functions as ff,
                       iterable_functions as itf,
                       sequence_functions as sf,
                       blast_functions as bf)
except ModuleNotFoundError:
    from SchemaRefinery.utils import (file_functions as ff,
                                      iterable_functions as itf,
                                      sequence_functions as sf,
                                      blast_functions as bf)

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
		locus_id = ff.file_basename(locus, file_extension=False)

		# Create paths to gene files in new schema
		locus_file = ff.join_paths(schema_path, [f'{locus_id}.fasta'])
		locus_short_file = ff.join_paths(schema_short_path, [f'{locus_id}_short.fasta'])

		# Create temp directory for current gene
		locus_temp_dir = ff.join_paths(schema_path, [f'{locus_id}_temp'])
		# Create temp directory for the current gene
		ff.create_directory(locus_temp_dir)

		# Dictionaries mapping gene identifiers to DNA and Protein sequences
		locus_seqs = sf.import_sequences(locus)
		prot_seqs = {k: str(sf.translate_sequence(v, table_id)) for k, v in locus_seqs.items()}

		if len(locus_seqs) > 1:
			# Identify DNA sequences that code for same protein
			equal_prots = sf.determine_duplicated_seqs(prot_seqs)

			# Get only one identifier per protein
			ids_to_blast = [protids[0] for protein, protids in equal_prots.items()]

			# Get longest sequence as first representative
			longest = sf.determine_longest(ids_to_blast, prot_seqs)
			representatives.append(longest)
			final_representatives.append(longest)

			# Create FASTA file with distinct protein sequences
			protein_file = ff.join_paths(locus_temp_dir, [f'{locus_id}_protein.fasta'])
			protein_data = [[i, prot_seqs[i]] for i in ids_to_blast]
			protein_lines = sf.fasta_lines('>{0}\n{1}', protein_data)
			ff.write_lines(protein_lines, protein_file)

			# Create blastdb with all distinct proteins
			blastp_db = ff.join_paths(locus_temp_dir, [locus_id])
			db_std = bf.make_blast_db(makeblastdb_path, protein_file, blastp_db, 'prot')

			# Determine appropriate blastp task (proteins < 30aa need blastp-short)
			blastp_task = bf.determine_blast_task(equal_prots)

			# Cycle BLAST representatives against non-representatives until
			# all non-representatives have a representative
			while len(set(ids_to_blast) - set(representatives)) != 0:
				# create FASTA file with representative sequences
				rep_file = ff.join_paths(locus_temp_dir, [f'{locus_id}_rep_protein.fasta'])
				rep_protein_data = [[r, prot_seqs[r]] for r in representatives]
				rep_protein_lines = sf.fasta_lines('>{0}\n{1}', rep_protein_data)
				ff.write_lines(rep_protein_lines, rep_file)

				# Compute self-score for representative alleles
				for seqid in representatives:
					if seqid not in rep_self_scores:
						record = sf.fasta_str_record('>{0}\n{1}', [seqid, prot_seqs[seqid]])
						current_rep_file = ff.join_paths(locus_temp_dir, [f'{seqid}_solo.fasta'])
						ff.write_lines([record], current_rep_file)
						# Create file with representative seqid to only compare against self
						id_file = ff.join_paths(locus_temp_dir, [f'{seqid}_ids.txt'])
						ff.write_lines([seqid], id_file)
						binary_file = f'{id_file}.bin'
						blast_std = bf.run_blastdb_aliastool(blastdb_aliastool_path, id_file, binary_file)
						id_file = binary_file

						rep_blastout = ff.join_paths(locus_temp_dir, [f'{seqid}_blastout.tsv'])
						# Cannot get self-alignemnt for some sequences if composition-based stats is enabled
						blast_std = bf.run_blast(blastp_path, blastp_db, current_rep_file,
											  rep_blastout, 1, 1,
											  id_file, 'blastp', None, 0)
						rep_results = ff.read_tabular(rep_blastout)
						if len(rep_results) > 0:
							rep_self_scores[rep_results[0][0]] = float(rep_results[0][6])
						else:
							print('Could not determine the self-alignment raw '
								f'score for {rep_results[0][0]}')

				# Create file with seqids to BLAST against
				ids_str = itf.join_list([str(i) for i in ids_to_blast if i not in representatives], '\n')
				ids_file = ff.join_paths(locus_temp_dir, [f'{locus_id}_ids.txt'])
				ff.write_to_file(ids_str, ids_file, 'w', '')
				binary_file = f'{ids_file}.bin'
				blast_std = bf.run_blastdb_aliastool(blastdb_aliastool_path, ids_file, binary_file)
				ids_file = binary_file

				# BLAST representatives against non-represented
				blast_output = ff.join_paths(locus_temp_dir, [f'{locus_id}_blast_out.tsv'])
				# Set 'max_target_seqs' to huge number because BLAST only
				# returns 500 hits by default
				blast_std = bf.run_blast(blastp_path, blastp_db, rep_file,
									  blast_output, 1, 1, ids_file,
									  blastp_task, 100000)

				# Import BLAST results
				blast_results = ff.read_tabular(blast_output)

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
				hotspot_reps = set(itf.flatten_list(list(hot_reps.values())))
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
		locus_lines = sf.fasta_lines('>{0}\n{1}', locus_data)
		ff.write_lines(locus_lines, locus_file)

		# Get total number of valid sequences
		total_sequences = len(locus_lines)

		# Write schema file with representatives
		locus_rep_data = [[r, locus_seqs[r]]
						  for r in final_representatives]
		locus_rep_lines = sf.fasta_lines('>{0}\n{1}', locus_rep_data)
		ff.write_lines(locus_rep_lines, locus_short_file)

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
	schema_short_path = ff.join_paths(schema_path, ['short'])

	# Create output directories
	ff.create_directory(schema_path)
	ff.create_directory(schema_short_path)

	# Import list of loci to adapt
	loci_list = ff.read_lines(input_file, strip=True)
	print(f'Number of loci to adapt: {len(loci_list)}')
	# Count number of sequences and mean length per locus
	loci_info = []
	loci_pools = multiprocessing.Pool(processes=cpu_cores)
	gp = loci_pools.map_async(sf.fasta_stats, loci_list, callback=loci_info.extend)
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
