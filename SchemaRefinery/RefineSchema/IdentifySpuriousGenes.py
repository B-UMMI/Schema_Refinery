#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import csv
import sys
import shutil
import subprocess
from typing import Any, List, Dict, Optional, Set, Tuple

import pandas as pd

try:
	from SchemaAnnotation import (consolidate as cs)
	from utils import (core_functions as cof,
						clustering_functions as cf,
										file_functions as ff,
										alignments_functions as af,
										sequence_functions as sf,
										iterable_functions as itf,
										blast_functions as bf,
										linux_functions as lf,
										classify_cds_functions as ccf,
										constants as ct,
										Types as tp,
										print_functions as pf,
										logger_functions as logf,
										time_functions as ti,
										globals as gb)
except ModuleNotFoundError:
	from SchemaRefinery.SchemaAnnotation import (consolidate as cs)
	from SchemaRefinery.utils import (core_functions as cof,
										clustering_functions as cf,
										file_functions as ff,
										alignments_functions as af,
										sequence_functions as sf,
										iterable_functions as itf,
										blast_functions as bf,
										linux_functions as lf,
										classify_cds_functions as ccf,
										constants as ct,
										Types as tp,
										print_functions as pf,
										logger_functions as logf,
										time_functions as ti,
										globals as gb)


def compute_cds_frequency(input_fasta, temp_directory):
	"""
	Calculates the CDS frequency in the genomes in 3 different ways depending on the run mode.
	"""
	# Import unclassified CDSs (sequence IDs as keys and sequences as values)
	pf.print_message("Importing unclassified CDSs...", "info")
	unclassified_cds = sf.fetch_fasta_dict(input_fasta)
	# Count the frequency of CDSs in the genomes and identify their presence
	pf.print_message("Identifying CDSs present in the schema and counting frequency of missing CDSs in the genomes...", "info")
	cds_present: str = os.path.join(temp_directory, "2_cds_preprocess/cds_deduplication/distinct.hashtable")
	decoded_seqids: Dict[str, List[float]] = itf.decode_CDS_sequences_ids(cds_present)
	frequencies = ccf.count_cds_frequency(unclassified_cds, decoded_seqids)

	return frequencies


def compute_loci_frequency(input_directory):
	"""
	"""
	# List FASTA files in the schema folder and map them to their full paths
	profiles_file = ff.join_paths(input_directory, ["results_alleles.tsv"])
	profiles = pd.read_csv(profiles_file, sep = '\t', dtype = object)
	loci_ids = profiles.columns[1:]
	frequencies: Dict[str, int] = {}
	for locus in loci_ids:
		locus_column = profiles[locus]
		# Convert LNF, ASM, and ALM classifications to 0 (LNF, ASM, ALM) and everything else to 1 (assume NIPH/NIPHEM and PLOT5/PLOT3 as presence)
		masked_column: pd.Series = locus_column.map(lambda x: 0 if str(x) in {'LNF', 'ASM', 'ALM'} else 1)
		# Sum all 1's to get locus frequency in the genomes
		frequencies[locus] = masked_column.sum()

	return frequencies


def identify_spurious_genes(schema_directory: List[str], output_directory: str, allelecall_directory: List[str],
							annotations: str, constants: List[Any], run_mode: str, cpu: int, no_cleanup: bool) -> None:
	"""
	Identify spurious genes in the given schema.

	Parameters
	----------
	schema_directory : List[str]
		Path to the schema directory.
	output_directory : str
		Path to the output directory.
	allelecall_directory : List[str]
		Path to the allele call directory.
	annotations : str
		Path to a TSV file containing loci annotations.
	constants : List[Any]
		List of constants used in the process.
	run_mode : str
		Mode of running the process.
	cpu : int
		Number of CPUs to use.
	no_cleanup : bool
		Flag to indicate whether to clean up temporary files.

	Returns
	-------
	None
	"""
	# Create base output directory
	ff.create_directory(output_directory)
	# Create initial processing output directory based on run mode
	processed_indata: str = os.path.join(output_directory, '1_Input_Data')
	ff.create_directory(processed_indata)
	# Create folder to store BLAST inputs
	blast_inputs = os.path.join(output_directory, "2_BLAST_Inputs")
	ff.create_directory(blast_inputs)
	# Set minimum number of genomes in which a CDS cluster should be present to be kept in the analysis
	# Set to 5 if number of genomes <= 20 else to 1% of the dataset size
	if not constants[2]:
		pf.print_message("Setting value for the minimum number of genomes a CDS cluster should be present in...", "info")
		ccf.set_minimum_genomes_threshold(allelecall_directory[0], constants)
		pf.print_message(f"Set genome presence threshold to {constants[2]}.", "info")

	# Get Path to the BLAST executables
	makeblastdb_exec: str = lf.get_tool_path('makeblastdb')
	blastn_exec: str = lf.get_tool_path('blastn')
	blastp_exec: str = lf.get_tool_path('blastp')
	blastdb_aliastool_exec: str = lf.get_tool_path('blastdb_aliastool')

	excluded = []
	if run_mode == 'unclassified_cds':
		# Copy file with the unclassified CDSs to the initial processing output directory
		pf.print_message(f"Copying FASTA file with input CDSs to {processed_indata}...", "info")
		unclassified_cds_src: str = os.path.join(allelecall_directory[0], "unclassified_sequences.fasta")
		unclassified_cds_dest = shutil.copy(unclassified_cds_src, processed_indata)
		# Define the path to the temp folder with the intermediate files from AlleleCall
		temp_folder: str = os.path.join(allelecall_directory[0], 'temp')
		# Calculate the frequency of the CDSs in the genomes
		cds_frequency = compute_cds_frequency(unclassified_cds_dest, temp_folder)

		# Exclude CDSs shorter than size threshold
		pf.print_message(f"Excluding CDSs shorter than {constants[5]} bp...", "info")
		unclassified_cds = sf.fetch_fasta_dict(unclassified_cds_dest)
		unclassified_cds, short = ccf.filter_by_size(unclassified_cds, constants[5])
		excluded.append(short)
		pf.print_message(f"{len(unclassified_cds)}/{len(cds_frequency)} CDSs have size greater or equal to {constants[5]} bp.", 'info')

		pf.print_message("Translating and deduplicating the CDSs...", "info")
		# Translate and deduplicate the CDSs
		# The protein_sequences contains only the distinct sequences
		protein_sequences, _ = sf.translate_seq_deduplicate(unclassified_cds, processed_indata, constants[5], constants[6], True)

		pf.print_message("Clustering the translated distinct CDSs...", "info")
		# Sort proteins by size before clustering
		protein_sequences = {k: v for k, v in sorted(protein_sequences.items(), key=lambda x: len(x[1]), reverse=True)}
		# Run minimizer-based clustering
		reps_groups: Dict[str, List[str]] = {}
		clusters: Dict[str, List[str]] = {}
		reps_sequences: Dict[str, str] = {}
		prot_len_dict: Dict[str, int] = {}
		clusters, reps_sequences, reps_groups, prot_len_dict = cf.minimizer_clustering(
			protein_sequences,
			clusters,
			reps_sequences,
			reps_groups,
			clustering_sim=constants[3],
			clustering_cov=constants[4],
			size_threshold = constants[8],
			grow = True)

		# Reformat clusters to contain only representative IDs as keys and a list of CDS IDs as values
		clusters = {repid: [m[0] for m in members] for repid, members in clusters.items()}

		# Determine the total number of clusters
		total_clusters: int = len(clusters)
		pf.print_message(f"Clustered {len(protein_sequences)} distinct proteins into {total_clusters} clusters.", "info")

		# Determine the number of singleton clusters
		singletons: int = len([cluster for cluster in clusters.values() if len(cluster) == 1])
		pf.print_message(f"{total_clusters - singletons} clusters contain more than one CDS.", "info")
		pf.print_message(f"{singletons} clusters are singletons.", "info")

		pf.print_message(f"Excluding clusters based on the frequency of the CDSs in the genomes (<{constants[2]})...", "info")
		# Sum the frequency of the CDSs in each cluster to get the frequency per cluster
		frequency_in_genomes: Dict[str, int] = {repid: sum(cds_frequency[seqid] for seqid in members) for repid, members in clusters.items()}
		# Drop clusters with a frequency below the given threshold and add the corresponding CDS IDs to the dropped dict
		excluded_seqids: List[str] = []
		excluded_clusters: List[str] = []
		for repid in clusters:
			if frequency_in_genomes[repid] < constants[2]:
				excluded_seqids.extend(clusters[repid])
				excluded_clusters.append(repid)
		# Exclude clusters
		for repid in excluded_clusters:
			del clusters[repid]

		# Add seqids excluded based on cluster frequency on genomes
		excluded.append(excluded_seqids)

		pf.print_message(f"Excluded {len(excluded_clusters)} clusters (a total of {len(excluded_seqids)} CDSs).")
		pf.print_message(f"{len(clusters)} clusters remain after filtering based on CDS frequency.", "info")

		# Create the files needed to run BLAST
		pf.print_message("Creating input files for BLAST...", "info")
		loci_files, concatenated_files = cof.prepare_cluster_blast_infiles(blast_inputs, clusters, unclassified_cds, protein_sequences, blastdb_aliastool_exec)

	if run_mode == 'schema':
		# Copy the schema to the folder with intermediate files
		pf.print_message("Copying input schema to temp folder...", "info")
		temp_schema = ff.join_paths(processed_indata, [ff.file_basename(schema_directory[0])])
		ff.copy_folder(schema_directory[0], temp_schema)
		# Calculate the frequency of the loci in the genomes
		frequency_in_genomes = compute_loci_frequency(allelecall_directory[0])
		pf.print_message("Creating input files for BLAST...", "info")
		loci_files, concatenated_files = cof.prepare_loci_blast_infiles(temp_schema, blast_inputs, constants[5], constants[6], blastdb_aliastool_exec)

	if run_mode == 'schema_vs_schema':
		# Copy schemas to temp folder
		ff.copy_folder(schema_directory[0], processed_indata)
		ff.copy_folder(schema_directory[1], processed_indata)

		pf.print_message(f'First schema: {schema_directory[0]}', 'info')
		pf.print_message(f'Second schema: {schema_directory[1]}', 'info')

		# Calculate the loci frequency bases on the allele calling results with the first schema
		first_schema_frequencies = compute_loci_frequency(allelecall_directory[0])
		# Second schema
		second_schema_frequencies = compute_loci_frequency(allelecall_directory[1])
		frequency_in_genomes = [first_schema_frequencies, second_schema_frequencies]

		### What's this???
		ff.merge_folders(schema_directory[0], schema_directory[1], schema_folder, cpu, bsr, translation_table)

		pf.print_message('Prepare loci files for Blast and count frequencies.', 'info')
		(all_nucleotide_sequences,
		master_file_path,
		trans_paths,
		to_blast_paths,
		all_alleles,
		group_reps_ids,
		group_alleles_ids,
		to_run_against,
		new_max_hits,
		seqid_file_dict) = cof.prepare_loci(schema_folder,
											constants,
											initial_processing_output)

	#####################################
	# Main part common to all run modes #
	#####################################

	# Create output folder to store BLASTp results
	blast_output: str = os.path.join(output_directory, '3_BLAST_Results')
	ff.create_directory(blast_output)

	# Create BLASTp database
	pf.print_message(f"Creating BLASTp database for {concatenated_files[1]}...", "info")
	blastp_db_directory: str = ff.join_paths(blast_inputs, ['BLASTp_db'])
	ff.create_directory(blastp_db_directory)
	blastp_db_path: str = ff.join_paths(blastp_db_directory, ['BLASTp_db'])
	bf.make_blast_db(makeblastdb_exec, concatenated_files[1], blastp_db_path, 'prot')

	# Calculate BLASTp self-score for all loci/cluster representatives
	self_blastp_inputs = {locus: paths[0] for locus, paths in loci_files.items()}
	self_blastp_outdir = ff.join_paths(blast_output, ["self_scores"])
	ff.create_directory(self_blastp_outdir)
	self_score_dict: Dict[str, float] = bf.calculate_self_score(self_blastp_inputs, blastp_exec, self_blastp_outdir, cpu)

	# Run BLASTp to align each representative protein against all other representative proteins
	pf.print_message("Running BLASTp...", "info")
	blastp_inputs: Dict[str, List[str]] = {locus: [paths[1], paths[3]] for locus, paths in loci_files.items()}
	blastp_outdir = ff.join_paths(blast_output, ["BLASTp_results"])
	ff.create_directory(blastp_outdir)
	# Allow multiple High-Scoring Segment Pairs to determine the global aligned fraction
	# Define a value for max_targets that is high enough to report all/most relevant alignments
	blastp_outfiles = bf.run_blast_operations(blastp_inputs, blastp_outdir, blastp_db_path, blastp_exec, cpu, max_hsps=5, max_targets=100)

	# Create BLASTn database
	blastn_db_directory: str = ff.join_paths(blast_inputs, ['BLASTn_db'])
	ff.create_directory(blastn_db_directory)
	blastn_db_path: str = ff.join_paths(blastn_db_directory, ['BLASTn_db'])
	bf.make_blast_db(makeblastdb_exec, concatenated_files[0], blastn_db_path, 'nucl')

	# Run BLASTn
	# Only using representative alleles, so running for all loci should not be a problem
	pf.print_message("Running BLASTn...", "info")
	blastn_inputs: Dict[str, List[str]] = {locus: [paths[0], paths[3]] for locus, paths in loci_files.items()}
	blastn_outdir = ff.join_paths(blast_output, ["BLASTn_results"])
	ff.create_directory(blastn_outdir)
	blastn_outfiles = bf.run_blast_operations(blastn_inputs, blastn_outdir, blastn_db_path, blastn_exec, cpu, max_hsps=5, max_targets=100)

	# Process BLAST results
	blast_results_dict: Dict[str, Dict[str, List[List[str], List[str]]]] = {}
	for i, file in enumerate(blastp_outfiles):
		blastp_file = file
		# Import BLASTp results
		with open(file, 'r') as infile:
			# Each line contains the following values: qseqid sseqid qlen slen qstart qend sstart send length score gaps pident
			blastp_lines = list(csv.reader(infile, delimiter='\t'))
			blastp_data = {}
			if blastp_lines:
				# Group all HSPSs for the same query-subject pair
				for line in blastp_lines:
					blastp_data.setdefault((line[0], line[1]), []).append(line[2:])

		# Import BLASTn results
		blastn_file = blastn_outfiles[i]
		with open(blastn_file, 'r') as infile:
			# Each line contains the following values: qseqid sseqid qlen slen qstart qend sstart send length score gaps pident
			blastn_lines = list(csv.reader(infile, delimiter='\t'))
			blastn_data = {}
			if blastn_lines:
				# Group all HSPSs for the same query-subject pair
				for line in blastn_lines:
					blastn_data.setdefault((line[0], line[1]), []).append(line[2:])

		# No matches with both BLASTp and BLASTn
		if not blastp_data and not blastn_data:
			pf.print_message(f"No BLASTp and BLASTn matches for {file}.", "info")
			continue

		for qs_pair, matches in blastp_data.items():
			query, subject = qs_pair
			# Determine and add global values
			if qs_pair not in blast_results_dict and (subject, query) not in blast_results_dict:
				# Get query and subject frequency values to the line
				# Use the seqid in the BLAST output files if running in the "unclassified_cds" mode
				# Determine the locus identifier otherwise
				query_frequency = frequency_in_genomes.get(query, frequency_in_genomes[query.rsplit('_', 1)[0]])
				subject_frequency = frequency_in_genomes.get(subject, frequency_in_genomes[subject.rsplit('_', 1)[0]])
				if run_mode == "schema_vs_schema":
					query_frequency_ss = frequency_in_genomes_second_schema.get(query, frequency_in_genomes[query.rsplit('_', 1)[0]])
					subject_frequency_ss = frequency_in_genomes_second_schema.get(subject, frequency_in_genomes[subject.rsplit('_', 1)[0]])
				else:
					query_frequency_ss = None
					subject_frequency_ss = None

				# Calculate the frequency ratio
				# The frequency of the query or subject should not be 0 because of the frequency threshold initially applied
				frequency_ratio = round(min(query_frequency/subject_frequency, subject_frequency/query_frequency), 3)

				# Get the sequence length for the query and subject
				query_len = matches[0][0]
				subject_len = matches[0][1]

				global_data = [query_len, subject_len, query_frequency, subject_frequency, query_frequency_ss, subject_frequency_ss, frequency_ratio]
				# Create query-subject pair entry
				blast_results_dict.setdefault(qs_pair, [[[],[]], [[],[]], global_data])

			# Get data for each match
			for m in matches:
				bsr = round(float(m[7]) / self_score_dict[query], 3)
				# Get the length of the alignment (have to merge the length of the HSPs if there are multiple HSPs for the same query-subject pair later)
				query_aligned_interval = [int(m[2]), int(m[3])] if int(m[2]) < int(m[3]) else [int(m[3]), int(m[2])]
				subject_aligned_interval = [int(m[4]), int(m[5])] if int(m[4]) < int(m[5]) else [int(m[5]), int(m[4])]
				match_pident = m[9]
				match_data = [bsr, query_aligned_interval, subject_aligned_interval, match_pident]
				if qs_pair in blast_results_dict:
					blast_results_dict[qs_pair][0][0].append(match_data)
				elif (subject, query) in blast_results_dict:
					blast_results_dict[(subject, query)][1][0].append(match_data)

			# Get BLASTn data
			if qs_pair in blastn_data or (subject, query) in blastn_data:
				blastn_matches = blastn_data.get(qs_pair, (subject, query))
				for m in blastn_matches:
					bsr = None
					query_aligned_interval = [int(m[2]), int(m[3])] if int(m[2]) < int(m[3]) else [int(m[3]), int(m[2])]
					subject_aligned_interval = [int(m[4]), int(m[5])] if int(m[4]) < int(m[5]) else [int(m[5]), int(m[4])]
					match_pident = m[9]
					match_data = [bsr, query_aligned_interval, subject_aligned_interval, match_pident]
					if qs_pair in blast_results_dict:
						blast_results_dict[qs_pair][0][1].append(match_data)
					elif (subject, query) in blast_results_dict:
						blast_results_dict[(subject, query)][1][1].append(match_data)

#########################

	# Compute the total alignment length for each query-subject pair by merging the HSPs if there are multiple HSPs for the same query-subject pair
	# At this stage, the results are still structured as Query-Subject pairs of representative alleles
	# It is necessary to merge the results by locus at a later stage
	for qs_pair, match_data in blast_results_dict.items():
		query_length: int = int(match_data[2][0])
		subject_length: int = int(match_data[2][1])
		# Merge aligned intervals
		# Query to subject alignments
		# Query
		query_aligned_intervals: List[int] = sorted([m[1] for m in match_data[0][0]], key= lambda x: x[0])
		query_aligned_merged_intervals = af.merge_intervals(query_aligned_intervals)
		query_aligned_length = sum([(i[1]-i[0]) for i in query_aligned_merged_intervals])
		# Subject
		subject_aligned_intervals: List[int] = sorted([m[2] for m in match_data[0][0]], key= lambda x: x[0])
		subject_aligned_merged_intervals = af.merge_intervals(subject_aligned_intervals)
		subject_aligned_length = sum([(i[1]-i[0]) for i in subject_aligned_merged_intervals])
		# Compute the weighted pident based on the pident of all HSPSs
		# Different HSPSs can overlap and it is not possible to compute the true pident based on BLAST's outmft 6
		query_aligned_identity: float = sum(float(m[3])*((m[1][1]-m[1][0])/query_length) for m in match_data[0][0])
		subject_aligned_identity: float = sum(float(m[3])*((m[2][1]-m[2][0])/subject_length) for m in match_data[0][0])
		# Compute the global palign minimum and maximum values
		global_talign_min: float = min(query_aligned_length / query_length,
									subject_aligned_length / subject_length)
		global_talign_max: float = max(query_aligned_length / query_length,
									subject_aligned_length / subject_length)
		# Compute the minimum and maximum global percent identity
		global_pident_min: float = min(query_aligned_identity, subject_aligned_identity)
		global_pident_max: float = max(query_aligned_identity, subject_aligned_identity)
		# Get the maximum BSR value from all query-subject HSPSs
		max_bsr: float = max([m[0] for m in match_data[0][0]])
		blast_results_dict[qs_pair][0].append([round(global_talign_min, 3), round(global_talign_max, 3),
											   round(global_pident_min, 3), round(global_pident_max, 3),
											   max_bsr])

		# Subject to query alignments
		# Query
		query_aligned_intervals: List[int] = sorted([m[1] for m in match_data[1][0]], key= lambda x: x[0])
		# It is possible that there are no Subject-Query alignments
		if len(query_aligned_intervals) > 0:
			query_aligned_merged_intervals = af.merge_intervals(query_aligned_intervals)
			query_aligned_length = sum([(i[1]-i[0]) for i in query_aligned_merged_intervals])
			# Subject
			subject_aligned_intervals: List[int] = sorted([m[2] for m in match_data[1][0]], key= lambda x: x[0])
			subject_aligned_merged_intervals = af.merge_intervals(subject_aligned_intervals)
			subject_aligned_length = sum([(i[1]-i[0]) for i in subject_aligned_merged_intervals])
			# Compute the weighted pident based on the pident of all HSPSs
			# Different HSPSs can overlap and it is not possible to compute the true pident based on BLAST's outmft 6
			query_aligned_identity: float = sum(float(m[3])*((m[1][1]-m[1][0])/query_length) for m in match_data[1][0])
			subject_aligned_identity: float = sum(float(m[3])*((m[2][1]-m[2][0])/subject_length) for m in match_data[1][0])
			# Compute the global palign minimum and maximum values
			global_talign_min: float = min(query_aligned_length / query_length,
										subject_aligned_length / subject_length)
			global_talign_max: float = max(query_aligned_length / query_length,
										subject_aligned_length / subject_length)
			# Compute the minimum and maximum global percent identity
			global_pident_min: float = min(query_aligned_identity, subject_aligned_identity)
			global_pident_max: float = max(query_aligned_identity, subject_aligned_identity)
			# Get the maximum BSR value from all query-subject HSPSs
			max_bsr: float = max([m[0] for m in match_data[1][0]])
			blast_results_dict[qs_pair][1].append([round(global_talign_min, 3), round(global_talign_max, 3),
												round(global_pident_min, 3), round(global_pident_max, 3),
												max_bsr])
		else:
			blast_results_dict[qs_pair][1].append([0, 0, 0, 0, 0])

#########################

	# Assign classes to all matches
	for qs_pair, match_data in blast_results_dict.items():
		query_freq = match_data[2][2]
		subject_freq = match_data[2][3]
		frequency_ratio = match_data[2][6]
		# Classify Query-Subject matches
		qs_data = match_data[0][2]
		# Classify based on global_talign_min and bsr
		if qs_data[0] >= constants[0]:
			if qs_data[4] >= constants[7]:
				assigned_class = '1a'
			else:
				assigned_class = '1b' if frequency_ratio <= 0.1 else '1c'
		elif (0.4 <= qs_data[0] < constants[0]):
			# Working with maximum pident
			if qs_data[3] >= constants[1]:
				if qs_data[1] >= constants[0]:
					assigned_class = '2a' if frequency_ratio <= 0.1 else '2b'
				else:
					assigned_class = '3a' if frequency_ratio <= 0.1 else '3b'
			else:
				if qs_data[1] >= constants[0]:
					assigned_class = '4a' if frequency_ratio <= 0.1 else '4b'
				else:
					assigned_class = '4c'
		else:
			assigned_class = '5'

		blast_results_dict[qs_pair][0].append(assigned_class)

		# Classify Subject-Query matches
		sq_data = match_data[1][2]
		# Classify based on global_talign_min and bsr
		if sq_data[0] >= constants[0]:
			if sq_data[4] >= constants[7]:
				assigned_class = '1a'
			else:
				assigned_class = '1b' if frequency_ratio <= 0.1 else '1c'
		elif (0.4 <= sq_data[0] < constants[0]):
			# Working with maximum pident
			if sq_data[3] >= constants[1]:
				if sq_data[1] >= constants[0]:
					assigned_class = '2a' if frequency_ratio <= 0.1 else '2b'
				else:
					assigned_class = '3a' if frequency_ratio <= 0.1 else '3b'
			else:
				if sq_data[1] >= constants[0]:
					assigned_class = '4a' if frequency_ratio <= 0.1 else '4b'
				else:
					assigned_class = '4c'
		else:
			assigned_class = '5'

		blast_results_dict[qs_pair][1].append(assigned_class)

	########## BLASTn matches are added as class 6 when there were no BLASTp matches and CDSs/loci without matches are classified as 7
	########## It's probably best if there's a match both with BLASTp and BLASTn and determine the difference to try to identify frameshifts (assume it's a frameshift if it only matches with BLASTn with high coverage)

#########################

	pf.print_message("Assigning actions based on classifications...", "info")
	assigned_classes = {}
	for qs_pair, match_data in blast_results_dict.items():
		query_class = match_data[0][-1]
		query_frequency = match_data[2][2]
		subject_class = match_data[1][-1]
		subject_frequency = match_data[2][3]
		# Need to group assigned classes for each locus if running in modes "schema" or "schema-vs-schema"
		if run_mode == "unclassified_cds":
			assigned_classes.setdefault(qs_pair[0], []).append([query_class, query_frequency, subject_frequency])
			assigned_classes.setdefault(qs_pair[1], []).append([subject_class, subject_frequency, query_frequency])
		else:
			query_id = (qs_pair[0]).rsplit('_', 1)[0]
			subject_id = (qs_pair[1]).rsplit('_', 1)[0]
			assigned_classes.setdefault(query_id, []).append([query_class, query_frequency, subject_frequency])
			assigned_classes.setdefault(subject_id, []).append([subject_class, subject_frequency, query_frequency])

	# Sort the results for each query by class
	class_order = ct.CLASSES_OUTCOMES
	# Assign an integer value to each class based on order
	order_index = {value: index for index, value in enumerate(class_order)}
	# Sort list of assigned classes
	actions = {}
	for seqid, classes in assigned_classes.items():
		sorted_classes = sorted(classes, key=lambda x: order_index.get(x[0]))
		top_class, query_frequency, subject_frequency = sorted_classes[0]
		# Assign action based on top class
		if top_class == '1a':
			actions[seqid] = ("Join", top_class)
		elif top_class in {'1b', '2a', '3a', '4a'}:
			# Need query and subject frequency to decide which to drop
			if query_frequency > subject_frequency:
				actions[seqid] = ("Add", top_class)
			elif subject_frequency > query_frequency:
				actions[seqid] = ("Drop", top_class)
		elif top_class in {'1c', '2b', '3b', '4b'}:
			actions[seqid] = ("Choice", top_class)
		# Classes 4c, 5, 6, and 7 are Add
		else:
			actions[seqid] = ("Add", top_class)

	# Use all Quey-Subject pairs to create a graph and identify connected components
	import networkx as nx
	G = nx.Graph()
	G.add_edges_from(blast_results_dict.keys())
	# Get connected components to define Join groups
	connected_components: List[Set[str]] = list(nx.connected_components(G))
	# Sort components by decreasing size
	connected_components = sorted(connected_components, key=lambda x: len(x), reverse=True)
	# Add action and top class to include in main output file
	groups = []
	for cc in connected_components:
		# Get the action for the unclassified CDS of for the locus if running the "schema" or "schema_vs_schema" mode
		if run_mode == "unclassified_cds":
			current_group = [(seqid, *actions[seqid]) for seqid in cc]
		else:
			seqids = set([seqid.rsplit('_', 1)[0] for seqid in cc])
			current_group = [(seqid, *actions[seqid]) for seqid in seqids]
		groups.append(current_group)

	# Add remaining representative seqids as Add
	if run_mode == "unclassified_cds":
		representative_seqids = clusters.keys()
		# Determine seqids that are not in any component/group
		component_seqids = set()
		for cc in connected_components:
			component_seqids.update(cc)
	# Add loci that had no matches
	else:
		representative_seqids = loci_files.keys()
		print(len(representative_seqids))
		component_seqids = set()
		for cc in connected_components:
			component_seqids.update(set([i.rsplit('_', 1)[0] for i in cc]))

	singletons = [seqid for seqid in representative_seqids if seqid not in component_seqids]
	pf.print_message(f"{len(singletons)} representative CDSs/loci had no matches with any sequence.", "info")
	for seqid in singletons:
		groups.append([(seqid, "Add", "7")])

	# Write main output file
	recommendations_file = ff.join_paths(output_directory, ["recommendations.tsv"])
	outlines = ["Locus\tAction\tClass"]
	i = 1
	for g in groups:
		# Add group number
		outlines.append(f"#{i}")
		lines = ["\t".join(a) for a in g]
		outlines.extend(lines)
		i += 1

	outtext = "\n".join(outlines)
	with open(recommendations_file, "w") as outfile:
		outfile.write(outtext+"\n")

	# Print the classification results
	# Count classifications
	action_counts = {}
	for g in groups:
		actions = [(a[1], a[2]) for a in g]
		for a in actions:
			if a in action_counts:
				action_counts[a] += 1
			else:
				action_counts[a] = 1

	sorted_counts = sorted(action_counts.items(), key=lambda x: order_index.get(x[0][1]))
	for scount in sorted_counts:
		pf.print_message(f"{scount[1]} representative sequences were classified as {scount[0][1]} and with an action of {scount[0][0]}.", "info")

	# Append loci annotations to the recommendations
	annotated_recommendations = os.path.join(output_directory, "recommendations_annotated.tsv")
	if annotations:
		pf.print_message(f"Appending annotations in {annotations} to file with recommendations...", "info")
		files = [recommendations_file, annotations]
		consolidated_annotations: str = cs.consolidate_annotations(files, False, annotated_recommendations)

	# Clean up temporary files
	if not no_cleanup:
		pf.print_message("Cleaning up temporary files...", "info")
		# Remove temporary files
		shutil.rmtree(processed_indata)
		shutil.rmtree(blast_inputs)
		shutil.rmtree(blast_output)


def main(schema_directory: List[str], output_directory: str, allelecall_directory: List[str],
		annotations: str, alignment_ratio_threshold: float, pident_threshold: float,
		clustering_sim_threshold: float, clustering_cov_threshold:float, genome_presence: int,
		absolute_size: int, translation_table: int, bsr: float, size_ratio: float, run_mode: str,
		cpu: int, no_cleanup: bool) -> None:
	"""
	Main function to identify spurious genes in a schema.

	Parameters
	----------
	schema_directory : List[str]
		Path to the schema directory.
	output_directory : str
		Path to the output directory.
	allelecall_directory : List[str]
		Path to the allele call directory.
	annotations : str
		Path to a TSV file containing loci annotations.
	alignment_ratio_threshold : float
		Threshold for alignment ratio.
	pident_threshold : float
		Threshold for percentage identity.
	clustering_sim_threshold : float
		Similarity threshold for clustering.
	clustering_cov_threshold : float
		Coverage threshold for clustering.
	genome_presence : int
		Minimum genome presence required.
	absolute_size : int
		Absolute size threshold.
	translation_table : int
		Genetic code used for translation.
	bsr : float
		BLAST Score Ratio value.
	size_ratio : float
		Size ratio threshold.
	run_mode : str
		Mode of running the process.
	cpu : int
		Number of CPU cores to use.
	no_cleanup : bool
		Flag to indicate whether to clean up temporary files.

	Returns
	-------
	None
	"""
	# Put all constants in one dict in order to decrease number of variables used around
	constants: List[Any] = [alignment_ratio_threshold, 
							pident_threshold,
							genome_presence,
							clustering_sim_threshold,
							clustering_cov_threshold,
							absolute_size,
							translation_table,
							bsr,
							size_ratio]

	identify_spurious_genes(schema_directory,
							output_directory,
							allelecall_directory,
							annotations,
							constants,
							run_mode,
							cpu,
							no_cleanup)
