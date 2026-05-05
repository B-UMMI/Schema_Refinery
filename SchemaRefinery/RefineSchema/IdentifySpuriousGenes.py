#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import csv
import sys
import subprocess
import pandas as pd
import shutil
from typing import Any, List, Dict, Optional, Set, Tuple


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


def compute_frequency(locus_column: pd.Series, locus_id: str) -> int:
	"""
	"""
	# For each locus count the frequency of ...
	# Convert LNF, ASM, and ALM classifications to 0 (LNF, ASM, ALM) and everything else to 1
	masked_column: pd.Series = locus_column.applymap(lambda x: 0 if str(x) == 'LNF' or str(x) == 'ASM' or str(x) == 'ALM' else 1)
	# Count the frequency of each locus in all genomes
	frequency = masked_column.sum()
	return frequency


def calculate_frequency(allelecall_directories: List[str], temp_paths: List[str], run_mode: str):
	"""
	Calculates the CDS frequency in the genomes in 3 different ways depending on the run mode.

	Parameters
	----------
	schema_folder : str
		Path to the folder containing schema FASTA files.
	allelecall_directory : str
		Path to the directory containing allele call results.
	second_schema_folder : str
		Path to the folder containing schema FASTA files of the second schema.
	second_allelecall_directory : str
		Path to the directory containing allele call results for the second schema.
	temp_paths : List[str]
		List of temporary paths.
	run_mode : str
		Mode for running the module.
		schema, unclassified_cds, schema_vs_schema

	Returns
	-------
	Tuple
		- frequency_in_genomes (Dict[str, int]): Dictionary of loci frequencies in genomes.
		- frequency_in_genomes_second_schema (Dict[str, int]): Dictionary of loci frequencies in genomes from the second given schema.
		- frequency_cds (Dict[str, int]): Dictionary of cds frequencies in genomes.
		- cds_presence_in_genomes (Dict[str, List[str]]): Dictionary of the presence of the cds in the genomes.
		- all_nucleotide_sequences (Dict[str, str]): Dictionary with CDS IDs as keys and sequences as values.
	"""
	# Calculate the CDS frequencies depending of the run mode
	pf.print_message(f"Calculating CDS frequency for the following run mode: {run_mode}", 'info')
	# Compute freuquencies for the schema and schema_vs_schema run modes
	if run_mode != "unclassified_cds":
		# List FASTA files in the schema folder and map them to their full paths
		fs_profiles_file = ff.join_paths(allelecall_directories[0], ["results_alleles.tsv"])
		fs_profiles = pd.read_csv(fs_profiles_file, sep = '\t', dtype = object)
		fs_loci_ids = fs_profiles.columns[1:]
		fs_frequencies: Dict[str, int] = {}
		for locus in fs_loci_ids:
			fs_frequencies[locus] = compute_frequency(fs_profiles[locus], locus)
		# Calculate frequencies for the second schema if run mode is schema_vs_schema
		if run_mode == "schema_vs_schema":
			ss_profiles_file = ff.join_paths(allelecall_directories[1], ["results_alleles.tsv"])
			ss_profiles = pd.read_csv(ss_profiles_file, sep = '\t', dtype = object)
			ss_loci_ids = ss_profiles.columns[1:]
			ss_frequencies: Dict[str, int] = {}
			for locus in ss_loci_ids:
				ss_frequencies[locus] = compute_frequency(ss_profiles[locus], locus)
		return fs_frequencies, ss_frequencies if run_mode == "schema_vs_schema" else fs_frequencies, None
	if run_mode == "unclassified_cds":
		# Import unclassified CDSs (sequence IDs as keys and sequences as values)
		pf.print_message("Importing unclassified CDSs...", "info")
		unclassified_cds = sf.fetch_fasta_dict(temp_paths[0])
		# Count the frequency of CDSs in the genomes and identify their presence
		pf.print_message("Identifying CDSs present in the schema and counting frequency of missing CDSs in the genomes...", "info")
		cds_present: str = os.path.join(temp_paths[1], "2_cds_preprocess/cds_deduplication/distinct.hashtable")
		decoded_seqids: Dict[str, List[float]] = itf.decode_CDS_sequences_ids(cds_present)
		cds_frequency = ccf.count_cds_frequency(unclassified_cds, decoded_seqids)
		return cds_frequency


# schema_directory = "/home/rmamede/WORK/TestSR/IdentifySpuriousGenes/mpneumoniae_wgMLST_schema_15042026"
# output_directory = "/home/rmamede/WORK/TestSR/IdentifySpuriousGenes/test_output"
# allelecall_directory = ["/home/rmamede/WORK/TestSR/IdentifySpuriousGenes/atb_plus_ncbi_results"]
# annotation_paths = []
# constants = [0.9, 90, None, 0.9, 0.9, 201, 4, 0.6, 0.8]
# run_mode = 'unclassified_cds'
# cpu = 4
# no_cleanup = True
def identify_spurious_genes(schema_directory: List[str], output_directory: str, allelecall_directory: List[str],
							annotation_paths: List[str], constants: List[Any], run_mode: str, cpu: int, no_cleanup: bool) -> None:
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
	annotation_paths : List[str]
		Paths for the files with the annotations.
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
	processed_indata: str = os.path.join(output_directory, '1_Initial_Data')
	ff.create_directory(processed_indata)
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

	if run_mode == 'unclassified_cds':
		# COpy file with the unclassified CDSs to the initial processing output directory
		unclassified_cds_src: str = os.path.join(allelecall_directory[0], "unclassified_sequences.fasta")
		unclassified_cds_dest = shutil.copy(unclassified_cds_src, processed_indata)
		# Deefine the path to the temp folder with the intermediate files from AlleleCall
		temp_folder: str = os.path.join(allelecall_directory[0], 'temp')
		# Calculate the frequency of the CDSs in the genomes
		cds_frequency = calculate_frequency([None, None], [unclassified_cds_dest, temp_folder], run_mode)

		pf.print_message("Excluding short CDSs...", "info")
		unclassified_cds = sf.fetch_fasta_dict(unclassified_cds_dest)
		unclassified_cds, shorter = ccf.filter_by_size(unclassified_cds, constants[5])

		# Save the unclassified CDSs
		pf.print_message("Saving unclassified CDSs to file...", "info")
		unclassified_cds_outfile: str = os.path.join(processed_indata, 'unclassified_cds_noShort.fasta')
		ccf.write_cds_to_fasta(unclassified_cds, unclassified_cds_outfile)
		pf.print_message(f"{len(unclassified_cds)}/{len(cds_frequency)} CDSs have size greater or equal to {constants[5]} bp.", 'info')

		pf.print_message("Translating and deduplicating unclassified CDSs...", "info")
		# Translate and deduplicate the CDSs
		# The protein_sequences contains only the distinct sequences and the protein_hashes contains the mapping of the distinct sequences to the original CDS IDs (including duplicates)
		protein_sequences, protein_hashes = sf.translate_seq_deduplicate(unclassified_cds, processed_indata, constants[5], constants[6], True)

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

		# Calculate the total number of clusters
		total_clusters: int = len(clusters)
		pf.print_message(f"{len(protein_sequences)} distinct proteins have been clustered into {total_clusters} clusters.", "info")

		# Calculate the number of singleton clusters
		singletons: int = len([cluster for cluster in clusters.values() if len(cluster) == 1])
		pf.print_message(f"Out of those clusters, {total_clusters - singletons} have more than one CDS.", "info")
		pf.print_message(f"Out of those clusters, {singletons} are singletons.", "info")

		pf.print_message("Excluding clusters based on CDS frequency...", "info")
		# Sum the frequency of the CDSs in each cluster to get the frequency per cluster
		frequency_in_genomes: Dict[str, int] = {repid: sum(cds_frequency[seqid] for seqid in members) for repid, members in clusters.items()}

		# Drop clusters that are present in less genomes than the given threshold and add the corresponding CDS IDs to the dropped dict with the reason for dropping them
		pf.print_message(f"Filtering clusters based on the frequency of the CDSs in the genomes (>= {constants[2]})...", "info")
		infrequent: List[str] = []
		excluded_clusters: List[str] = []
		for repid in clusters:
			if frequency_in_genomes[repid] < constants[2]:
				infrequent.extend(clusters[repid])
				excluded_clusters.append(repid)

		for repid in excluded_clusters:
			del clusters[repid]

		pf.print_message(f"Excluded {len(excluded_clusters)} clusters (a total of {len(infrequent)} CDSs).")
		pf.print_message(f"{len(clusters)} clusters remain after filtering based on CDS frequency.", "info")

		# Calculate k-mer similarity between the cluster representatives
		# pf.print_message("Computing k-mer similarity and coverage between cluster representatives...", "info")
		# representative_proteins = {repid: protein_sequences[repid] for repid in clusters.keys()}
		# representative_kmers_sim: Dict[str, float] = ccf.calculate_kmers_similarity(representative_proteins, reps_groups, prot_len_dict)

		blast_inputs = os.path.join(output_directory, "BLAST_inputs")
		ff.create_directory(blast_inputs)
		# Create the files needed to run BLAST
		pf.print_message("Creating BLAST input files...", "info")
		loci_files, concatenated_files = ccf.prepare_files_to_blast(blast_inputs, clusters, unclassified_cds, protein_sequences, blastdb_aliastool_exec)

	if run_mode == 'schema':
		# Calculate the frequency of the cds in genomes
		results = calculate_frequency([allelecall_directory[0], None], None, None, None, run_mode, None)
		frequency_in_genomes, _, _, _, _, dropped_alleles, _ = results

		pf.print_message('Prepare loci files for Blast and count frequencies.', 'info')
		(all_nucleotide_sequences,
		master_file_path,
		trans_paths,
		to_blast_paths,
		all_alleles,
		group_reps_ids,
		group_alleles_ids,
		to_run_against,
		seqid_file_dict) = cof.prepare_loci(schema_directory[0], constants, initial_processing_output)

	if run_mode == 'schema_vs_schema':
		ff.copy_folder(schema_directory[0], schema_folder)
		second_schema_folder = os.path.join(initial_processing_output, 'second_schema')
		ff.create_directory(second_schema_folder)
		ff.copy_folder(schema_directory[1], second_schema_folder)

		pf.print_message(f'First schema: {schema_directory[0]}', 'info')
		pf.print_message(f'Second schema: {schema_directory[1]}', 'info')

		# Calculate the frequency of the cds in genomes for the first schema and second schema
		results = calculate_frequency([allelecall_directory[0], allelecall_directory[1]], None, None, None, run_mode, None)
		frequency_in_genomes, frequency_in_genomes_second_schema, _, _, _, dropped_alleles, _ = results

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

	# Main part that is the same for all run modes

	# Create output folder to store BLASTp results
	blast_output: str = os.path.join(output_directory, 'BLAST_results')
	ff.create_directory(blast_output)

	# Create BLASTp database
	pf.print_message("Creating BLASTp database...", "info")
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
	blastp_results_files = bf.run_blast_operations(blastp_inputs, blastp_outdir, blastp_db_path, blastp_exec, cpu)

####################################
	# Run BLASTn
	#### Only run for the loci/clusters without BLASTp matches?
	#### Maybe run this after BLASTp and include BLASTn results in the data structure created after processing the BLASTp results
	blastn_db_directory: str = ff.join_paths(blast_inputs, ['BLASTn_db'])
	ff.create_directory(blastn_db_directory)
	blastn_db_path: str = ff.join_paths(blastn_db_directory, ['BLASTn_db'])
	bf.make_blast_db(makeblastdb_exec, concatenated_files[0], blastn_db_path, 'nucl')

	# Run BLASTn
	blastn_inputs: Dict[str, List[str]] = {locus: [paths[0], paths[3]] for locus, paths in loci_files.items()}
	blastn_outdir = ff.join_paths(blast_output, ["BLASTn_results"])
	ff.create_directory(blastn_outdir)
	pf.print_message("Running BLASTn...", "info")
	blastn_results_files = bf.run_blast_operations(blastn_inputs, blastn_outdir, blastn_db_path, blastn_exec, cpu)

	# Process the BLASTn results
	# blastn_results_dict: Dict[str, Dict[str, List[str]]] = {}
	# for file in blastn_results_files:
	# 	with open(file, 'r') as infile:
	# 		#### Need to check if the number of results is equal to the number of max hits to adjust the max hits to N*N, where N is the number of sequences in the db?
	# 		# Each line contains the following values: qseqid sseqid qlen slen qstart qend sstart send length score gaps pident
	# 		file_lines = list(csv.reader(infile, delimiter='\t'))
	# 		for line in file_lines:
	# 			if line[0] not in blastn_results_dict:
	# 				blastn_results_dict[line[0]] = {line[1]: line}
	# 			else:
	# 				blastn_results_dict[line[0]][line[1]] = line

##########################
	print()

	# Process the BLASTp results
	####### Need to take into account that the results may include query-to-subject and subject-to-query alignments!
	blastp_results_dict: Dict[str, Dict[str, List[List[str], List[str]]]] = {}
	for file in blastp_results_files:
		with open(file, 'r') as infile:
			#### Need to check if the number of results is equal to the number of max hits to adjust the max hits to N*N, where N is the number of sequences in the db?
			# Each line contains the following values: qseqid sseqid qlen slen qstart qend sstart send length score gaps pident
			file_lines = list(csv.reader(infile, delimiter='\t'))

		# Do not continue if file is empty
		if not file_lines:
			continue

		queries = set()
		for line in file_lines:
			# Compute the BSR
			query = line[0]
			queries.add(query)
			subject = line[1]
			bsr = round(float(line[9]) / self_score_dict[query], 3)
			# Get the length of the alignment (have to merge the length of the HSPs if there are multiple HSPs for the same query-subject pair later)
			query_aligned_interval = [int(line[4]), int(line[5])] if int(line[4]) < int(line[5]) else [int(line[5]), int(line[4])]
			subject_aligned_interval = [int(line[6]), int(line[7])] if int(line[6]) < int(line[7]) else [int(line[7]), int(line[6])]
			# Add query and subject frequency values to the line
			query_frequency = frequency_in_genomes[query]
			subject_frequency = frequency_in_genomes[subject]
			if run_mode == 'schema_vs_schema':
				query_frequency_ss = frequency_in_genomes_second_schema[query]
				subject_frequency_ss = frequency_in_genomes_second_schema[subject]
			else:
				query_frequency_ss = None
				subject_frequency_ss = None

			match_data = line + [bsr, query_aligned_interval, subject_aligned_interval]
			global_data = [query_frequency, subject_frequency, query_frequency_ss, subject_frequency_ss]
			if query not in blastp_results_dict:
				blastp_results_dict[query] = {subject: [[match_data], *global_data]}
			else:
				if subject not in blastp_results_dict[query]:
					blastp_results_dict[query][subject] = [[match_data], *global_data]
				else:
					blastp_results_dict[query][subject][0].append(match_data)

		queries = list(queries)

		# Compute the total alignment length for each query-subject pair by merging the HSPs if there are multiple HSPs for the same query-subject pair
		##### Do this to compute all values for each file in a single pass and keep only the values necessary to assign the classes in the future
		#### What if there are multiple queries???
		for subject, match_data in blastp_results_dict[query].items():
			query_length: int = int(match_data[0][0][2])
			subject_length: int = int(match_data[0][0][3])
			# Merge aligned intervals
			# Query
			query_aligned_intervals: List[int] = sorted([match[13] for match in match_data[0]], key= lambda x: x[1])
			query_aligned_merged_intervals = af.merge_intervals(query_aligned_intervals)
			query_aligned_length = sum([(i[1]-i[0]) for i in query_aligned_merged_intervals])
			# Subject
			subject_aligned_intervals: List[int] = sorted([match[14] for match in match_data[0]], key= lambda x: x[1])
			subject_aligned_merged_intervals = af.merge_intervals(subject_aligned_intervals)
			subject_aligned_length = sum([(i[1]-i[0]) for i in subject_aligned_merged_intervals])

			# Compute the global percent identity for the query and subject
			#### Cannot really compute it like this, right?
			query_aligned_identity: float = sum(float(match[11])*((match[13][1]-match[13][0])/query_length) for match in match_data[0])
			subject_aligned_identity: float = sum(float(match[11])*((match[14][1]-match[14][0])/subject_length) for match in match_data[0])

			# Compute the global palign minimum and maximum values
			global_talign_min: float = min(query_aligned_length / query_length,
										subject_aligned_length / subject_length)
			global_talign_max: float = max(query_aligned_length / query_length,
										subject_aligned_length / subject_length)
			# Compute the minimum and maximum global percent identity
			global_pident_min: float = min(query_aligned_identity, subject_aligned_identity)
			global_pident_max: float = max(query_aligned_identity, subject_aligned_identity)

			# Compute the minimum local palign
			query_local_taligns: float = min([(i[1]-i[0])/query_length for i in query_aligned_intervals])
			subject_local_taligns: float = min([(i[1]-i[0])/subject_length for i in subject_aligned_intervals])
			local_talign_min: float = min(query_local_taligns, subject_local_taligns)

			# Assign class to query-subject pair
			# Calculate the frequency ratio
			query_freq: int = match_data[1]
			subject_freq: int = match_data[2]
			# Get the maximum BSR value from all query-subject matches
			max_bsr: float = max([match[12] for match in match_data[0]])

			# The query_freq or subjet_freq should not be 0 because of the frequency threshold initially applied
			freq_ratio = min(query_freq/subject_freq, subject_freq/query_freq)
			# Classify based on global_talign_min and bsr
			if global_talign_min >= constants[0]:
				if max_bsr >= constants[7]:
					assigned_class = '1a'
				elif freq_ratio <= 0.1:
					assigned_class = '1b' if freq_ratio <= 0.1 else '1c'
			elif 0.4 <= global_talign_min < constants[0]:
				#### Work with minimum or maximum pident???
				if global_pident_min >= constants[1]:
					if global_talign_max >= constants[0]:
						assigned_class = '2a' if freq_ratio <= 0.1 else '2b'
					else:
						assigned_class = '3a' if freq_ratio <= 0.1 else '3b'
				else:
					if global_talign_max >= constants[0]:
						assigned_class = '4a' if freq_ratio <= 0.1 else '4b'
					else:
						assigned_class = '4c'
			else:
				assigned_class = '5'

		########## BLASTn matches are added as class 6 when there were no BLASTp matches
		########## It's probably best if there's a match both with BLASTp and BLASTn and determine the difference to try to identify frameshifts (assume it's a frameshift if it only matches with BLASTn with high BSR)

			# Append computed values to the matches
			blastp_results_dict[query][subject].extend([query_aligned_length, subject_aligned_length,
														round(global_talign_min, 3), round(global_talign_max, 3),
														round(global_pident_min, 3), round(global_pident_max, 3),
														round(local_talign_min, 3), assigned_class])

		# Sort the results for each query by class
		class_order = ct.CLASSES_OUTCOMES
		# Assign an integer value to each class based on order
		order_index = {value: index for index, value in enumerate(class_order)}
		# Sort subject matches based on assigned class and add sorted data as new value to original dictionary
		sorted_matches = dict(sorted(blastp_results_dict[query].items(), key=lambda x: order_index.get(x[1][-1])))
		blastp_results_dict[query] = sorted_matches

##################################

	pf.print_message("Processing results to merge results and assign actions...", "info")
	##### Need to get the class count per query
	##### Assign a class to each match between a query and a subject or retain only the query-subject assignment?
	processed_results = {"Join": [set(), set()], "Choice": {}, "Drop": {}, "Add": {}}
	for query, subject_matches in blastp_results_dict.items():
		for subject, match_data in subject_matches.items():
			assigned_class = match_data[-1]
			query_freq = match_data[1]
			subject_freq = match_data[2]
			if assigned_class == '1a':
				# Add node data
				processed_results["Join"][0].update([query, subject])
				# Add edge data
				# Do not add inverse edges
				if (subject, query) not in processed_results["Join"][1]:
					processed_results["Join"][1].add((query, subject))
			elif assigned_class in {'1b', '2a', '3a', '4a'}:
				if query_freq > subject_freq:
					processed_results["Drop"].setdefault(subject, []).append(assigned_class)
				else:
					processed_results["Drop"].setdefault(query, []).append(assigned_class)
			elif assigned_class in {'1c', '2b', '3b', '4b'}:
				processed_results["Choice"].setdefault(query, []).append((subject, assigned_class))
			#### Everything else is Add?
			else:
				processed_results["Add"].setdefault(query, []).append(assigned_class)

	# Create graph from 1a matches and find connected components to create all Join groups
	# Only if there are any Join actions
	actions = {}
	if len(processed_results["Join"][1]) > 0:
		import networkx as nx
		G = nx.Graph()
		G.add_edges_from(processed_results["Join"][1])
		# Get connected components to define Join groups
		connected_components: List[Set[str]] = list(nx.connected_components(G))

		# Sort components by size to get larger components first

		# Start creating groups from components
		for c in components:
			actions[i] = []
			for seqid in c:
				# Change action from Join to Drop if seqid was assigned a Drop action
				if seqid in processed_results["Drop"]:
					actions[i].append([seqid, "Drop", processed_results["Drop"][seqid]])
					del processed_results["Drop"][seqid]
					continue
				# Remove from Choice actions if seqid was assigned a Join action
				if seqid in processed_results["Choice"]:
					del processed_results["Choice"][seqid]
				actions[i].append([seqid, "Join", "1a"])
			# Sort to get Join actions first
			actions[i] = sorted(actions[i], key= lambda x: x[1], reverse=True)
			i += 1

	# print(processed_results["Choice"])

	# Add remaining Choice actions
	for seqid, choices in processed_results["Choice"].items():
		# Add main query
		if seqid not in processed_results["Drop"]:
			actions[seqid] = [[seqid, "Choice", [choices[0][1]]]]
		else:
			actions[seqid] = [[seqid, "Drop", [choices[0][1]]]]
			del processed_results["Drop"][seqid]
		# Add matches
		for c in choices:
			if c[0] not in processed_results["Drop"]:
				actions[seqid].append([c[0], "Choice", [c[1]]])
			else:
				actions[seqid].append([c[0], "Drop", [c[1]]])
				del processed_results["Drop"][c[0]]

	# Add remaining Drop actions
	for seqid, assigned_class in processed_results["Drop"].items():
		actions[seqid] = [[seqid, "Drop", assigned_class]]

	# Everything else is added with the Add action
	for seqid, assigned_class in processed_results["Add"].items():
		if seqid not in processed_results["Choice"] and seqid not in processed_results["Drop"]:
			actions[seqid] = [[seqid, "Add", assigned_class]]

	print(actions)

	# Write main output file
	recommendations_file = ff.join_paths(output_directory, ["final_recommendations.tsv"])
	outlines = []
	for seqid, action_data in actions.items():
		for a in action_data:
			line = a[:2] + [a[2][0]]
			outlines.append("\t".join(line))
		outlines.append("#")
	
	outtext = "\n".join(outlines)
	with open(recommendations_file, "a") as outfile:
		outfile.write(outtext+"\n")

###############################################

	# count_results_by_class = itf.sort_subdict_by_tuple(count_results_by_class, classes_outcome)
	# # Extract which clusters are to maintain and to display to user.
	# clusters_to_keep: tp.MergedAllClasses
	# dropped_loci_ids: Set[str]
	# clusters_to_keep_1a, clusters_to_keep, dropped_loci_ids = cof.extract_clusters_to_keep(classes_outcome, count_results_by_class, drop_mark)

	# # Add the loci/new_loci IDs of the 1a joined clusters to the clusters_to_keep
	# clusters_to_keep_1a_renamed: Dict[str, List[str]] = {values[0]: values for key, values in clusters_to_keep_1a.items()}

	# # Merge classes
	# merged_all_classes: tp.MergedAllClasses = {'1a': clusters_to_keep_1a_renamed.copy()}
	# merged_all_classes.update(clusters_to_keep)
	# if run_mode == 'unclassified_cds':
	# 	updated_frequency_in_genomes: Dict[str, int] = ccf.update_frequencies_in_genomes(clusters_to_keep_1a_renamed,  frequency_in_genomes)
	# 	# Open dict to store IDs of the reps and alleles
	# 	group_reps_ids = {}
	# 	group_alleles_ids = {}
	# 	# Count the number of reps and alleles again because clusters were joined
	# 	group_reps_ids, group_alleles_ids = cof.count_number_of_reps_and_alleles(merged_all_classes,
	# 																			all_alleles,
	# 																			dropped_loci_ids,
	# 																			group_reps_ids,
	# 																			group_alleles_ids)

	# 	pf.print_message("Adding remaining clusters that didn't match by BLASTn...", "info")
	# 	# Add cluster not matched by BLASTn
	# 	all_matched_clusters: List[str] = itf.flatten_list([v for v in {key: value for key, value in clusters_to_keep.items() if key != '1a'}.values()]) + itf.flatten_list([values for values in clusters_to_keep_1a_renamed.values()])
	# 	clusters_to_keep['Retained_not_matched_by_blastn'] = set([cluster for cluster in all_alleles.keys() if cluster not in all_matched_clusters])

	# processed_drop: List[str] = []
	# # Add Ids of the dropped cases due to frequency during classification
	# cof.add_cds_to_dropped_cds(dropped_loci_ids,
	# 						dropped_alleles,
	# 						clusters_to_keep,
	# 						clusters_to_keep_1a_renamed,
	# 						all_alleles,
	# 						'Dropped_due_to_smaller_genome_presence_than_matched_cluster',
	# 						processed_drop)

	# pf.print_message("Extracting results...", "info")
	# related_clusters: tp.RelatedClusters
	# recommendations: tp.Recomendations
	# low_freq_recommendations: tp.Recomendations
	# # Extract the results from the processed results
	# related_clusters, recommendations, low_freq_recommendations = cof.extract_results(processed_results,
	# 																					count_results_by_class,
	# 																					frequency_in_genomes,
	# 																					frequency_in_genomes_second_schema if run_mode == 'schema_vs_schema' else None, 
	# 																					merged_all_classes,
	# 																					dropped_loci_ids,
	# 																					classes_outcome)

	# pf.print_message("Writing count_results_by_cluster.tsv, related_matches.tsv files and recommendations.tsv...", "info")
	# # Write the results to files and return the paths to the files.
	# reverse_matches: bool = True
	# (related_matches_path,
	#  count_results_by_cluster_path,
	#  recommendations_file_path,
	#  count_classes_final) = cof.write_recommendations_summary_results(to_blast_paths,
	# 																		related_clusters,
	# 																		count_results_by_class_with_inverse,
	# 																		group_reps_ids,
	# 																		group_alleles_ids,
	# 																		frequency_in_genomes,
	# 																		frequency_in_genomes_second_schema if run_mode == 'schema_vs_schema' else None,
	# 																		recommendations,
	# 																		low_freq_recommendations,
	# 																		reverse_matches,
	# 																		classes_outcome,
	# 																		output_directory)

	# # Get all of the CDS that matched with loci or loci matched with loci
	# is_matched: Dict[str, Any]
	# is_matched_alleles: Dict[str, Any]
	# is_matched, is_matched_alleles = cof.get_matches(all_relationships,
	# 												merged_all_classes,
	# 												sorted_blast_dict)

	# pf.print_message("Writing classes and cluster results to files...", "info")
	# report_file_path: str = os.path.join(blast_results, 'blast_all_matches.tsv')
	# # Write all of the alignments results to a file.
	# cof.alignment_dict_to_file(representative_blast_results,
	# 						   report_file_path,
	# 						   'w')
	# cof.alignment_dict_to_file(representative_blastn_results,
	# 						   report_file_path,
	# 						   'a')
	# # Write the processed results to a file alignments by clusters and classes.
	# cof.write_processed_results_to_file(merged_all_classes,
	# 								representative_blast_results,
	# 								classes_outcome,
	# 								all_alleles,
	# 								is_matched,
	# 								is_matched_alleles,
	# 								blast_results)
	# merged_classes_6 = {k: merged_all_classes[k] for k in {'6'}}
	# cof.write_processed_results_to_file(merged_classes_6,
	# 								representative_blastn_results,
	# 								['6'],
	# 								all_alleles,
	# 								is_matched,
	# 								is_matched_alleles,
	# 								blast_results)

	# if run_mode == 'unclassified_cds':
	# 	pf.print_message("Updating IDs and saving changes in cds_id_changes.tsv...", "info")
	# 	# Update the IDs and save the changes in a file.
	# 	ccf.update_ids_and_save_changes(merged_all_classes,
	# 								all_alleles,
	# 								cds_original_ids,
	# 								dropped_alleles,
	# 								all_nucleotide_sequences,
	# 								results_output)
	# 	pf.print_message("Writing dropped CDSs to file...", "info")
	# 	# Write the dropped CDS to a file.
	# 	ccf.write_dropped_cds_to_file(dropped_alleles, results_output)

	# pf.print_message("Writing dropped possible new loci to file...", "info")
	# # Write the dropped possible new loci to a file.
	# drop_possible_loci_output = cof.write_dropped_possible_new_loci_to_file(dropped_loci_ids,
	# 																	dropped_alleles,
	# 																	output_directory)
	# # Write alot of alleles file
	# alot_of_alleles_file = os.path.join(output_directory, 'alot_of_alleles.txt')
	# with open(alot_of_alleles_file, 'w') as alot_file:
	# 	for loci in loci_too_big:
	# 		alot_file.write(f"{loci}\n")

	# # Print the classification results
	# cof.print_classifications_results(merged_all_classes,
	# 									dropped_loci_ids,
	# 									to_blast_paths,
	# 									all_alleles,
	# 									count_classes_final)

	# # Annotate the recommendation outputs with the given annotation files using consolidate
	# pf.print_message("")
	# consolidated_annotations = os.path.join(output_directory, "recommendations_annotations.tsv")
	# if annotation_paths:
	# 	files: List[str]
	# 	files = [recommendations_file_path] + annotation_paths
	# 	pf.print_message("Consolidating annoations...", "info")
	# 	consolidated_annotations: str = cs.consolidate_annotations(files,
	# 								False,
	# 								consolidated_annotations)
	# 	pf.print_message('Annotation consolidation successfully completed.', 'info')
	# 	pf.print_message('')

	# # Clean up temporary files
	# if not no_cleanup:
	# 	pf.print_message("Cleaning up temporary files...", "info")
	# 	# Remove temporary files
	# 	ff.cleanup(output_directory, [related_matches_path,
	# 							alot_of_alleles_file,
	# 							count_results_by_cluster_path,
	# 							recommendations_file_path,
	# 							consolidated_annotations if annotation_paths else None,
	# 							drop_possible_loci_output,
	# 							temp_fastas_folder if run_mode == 'unclassified_cds' else None,
	# 							logf.get_log_file_path(gb.LOGGER)])


def main(schema_directory: List[str], output_directory: str, allelecall_directory: List[str],
		annotation_paths: List[str], alignment_ratio_threshold: float, pident_threshold: float,
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
	annotation_paths : List[str]
		Paths for the files with the annotations.
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
							annotation_paths,
							constants,
							run_mode,
							cpu,
							no_cleanup)
