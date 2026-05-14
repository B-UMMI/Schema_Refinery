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
	Calculate the CDS frequency in the genomes.

	Parameters
	----------
	input_fasta : str
		Path to the input FASTA file containing CDSs.
	temp_directory : str
		Path to the temp folder created by chewBBACA's AlleleCall
		module when the `--no-cleanup` option is provided.
 	"""
	# Import CDSs (sequence IDs as keys and sequences as values)
	cds_sequences = sf.fetch_fasta_dict(input_fasta)
	# Count the frequency of CDSs in the genomes and identify their presence
	# Read the file with the CDS genome presence data created by chewBBACA
	cds_presence: str = os.path.join(temp_directory, "2_cds_preprocess/cds_deduplication/distinct.hashtable")
	# Decode the data
	decoded_seqids: Dict[str, List[float]] = itf.decode_CDS_sequences_ids(cds_presence)
	# The CDSs and the decoded sequence IDs need to come from the same run
    # Each CDS must have a corresponding hash in the decoded sequence IDs
	frequencies: Dict[str, int] = {}
    for seqid, sequence in cds_sequences.items():
        hashed_seq: str = sf.seq_to_hash(str(sequence))
        # Count only the unique genome IDs for the frequency [1:]
        frequencies[seqid] = len(set(decoded_seqids[hashed_seq][1:]))

	return frequencies


def compute_class_counts(input_directory):
	"""
	"""
	profiles_file = ff.join_paths(input_directory, ["results_alleles.tsv"])
	profiles = pd.read_csv(profiles_file, sep = '\t', dtype = object)
	total_genomes = len(profiles["FILE"])
	loci_ids = profiles.columns[1:]
	classification_percentages = {}
	for locus in loci_ids:
		locus_column = profiles[locus]
		# Remove "INF-" prefixes
		masked_column: pd.Series = locus_column.str.replace('INF-', '')
		# Convert exact matches to 1
		masked_column: pd.Series = locus_column.map(lambda x: "1" if str(x) not in ct.CHEWIE_CLASSES else x)
		class_counts = masked_column.value_counts()
		# Determine percentage of special classifications
		class_counts = [class_counts.get(c, 0) for c in ct.CHEWIE_CLASSES] + [class_counts.get("1", 0)]
		class_counts = list(map(int, class_counts))
		# Do not include LNF classifications and exact matches
		special_percentage = round(sum(class_counts[:-2])/total_genomes, 3)
		# Determine frequency
		# Consider PLOT5/PLOT3/LOTSC and NIPH/NIPHEM classes as presence
		locus_frequency = sum(class_counts[2:7]+[class_counts[-1]])
		classification_percentages[locus] = [class_counts, special_percentage, locus_frequency]

	return classification_percentages


def identify_spurious_genes(input_schemas: List[str], output_directory: str, allelecall_directory: List[str],
							annotations: str, constants: List[Any], run_mode: str, cpu: int, no_cleanup: bool) -> None:
	"""
	Identify spurious genes in the given schema.

	Parameters
	----------
	input_schemas : List[str]
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
	# Set the minimum number of genomes in which a CDS cluster should be present to be kept in the analysis
	# Set to >=2% or 0 if dataset has less than 26 genomes
	if not constants[2]:
		pf.print_message("Setting value for the minimum number of genomes a representative CDS/locus should be present in...", "info")
		# The allele call results must be for the same dataset
		constants[2] = ccf.set_minimum_genomes_threshold(allelecall_directory[0], constants[2])
		pf.print_message(f"Set genome presence threshold to {constants[2]}.", "info")

	# Get Path to the BLAST executables
	makeblastdb_exec: str = lf.get_tool_path('makeblastdb')
	blastn_exec: str = lf.get_tool_path('blastn')
	blastp_exec: str = lf.get_tool_path('blastp')
	blastdb_aliastool_exec: str = lf.get_tool_path('blastdb_aliastool')

	excluded = {}
	if run_mode == 'cds':
		# Copy file with the unclassified CDSs to the initial processing output directory
		pf.print_message(f"Copying FASTA file with input CDSs to {processed_indata}...", "info")
		cds_sequences_src: str = os.path.join(allelecall_directory[0], "unclassified_sequences.fasta")
		cds_sequences_dest = shutil.copy(cds_sequences_src, processed_indata)
		# Define the path to the temp folder with the intermediate files from AlleleCall
		temp_folder: str = os.path.join(allelecall_directory[0], 'temp')
		# Calculate the frequency of the CDSs in the genomes
		pf.print_message("Computing the frequency of the CDSs in the genomes...", "info")
		cds_frequency = compute_cds_frequency(cds_sequences_dest, temp_folder)

		# Exclude CDSs shorter than size threshold
		pf.print_message(f"Excluding CDSs shorter than {constants[6]} bp...", "info")
		cds_sequences = sf.fetch_fasta_dict(cds_sequences_dest)
		cds_sequences, short = ccf.filter_by_size(cds_sequences, constants[6])
		for seqid in short:
			excluded.setdefault(seqid, []).append("short")

		pf.print_message(f"{len(cds_sequences)}/{len(cds_frequency)} CDSs have size greater or equal to {constants[6]} bp.", 'info')

		pf.print_message("Translating and deduplicating the CDSs...", "info")
		# Translate and deduplicate the CDSs
		# The protein_sequences contains only the distinct sequences
		protein_sequences, _ = sf.translate_seq_deduplicate(cds_sequences, processed_indata, constants[6], constants[7], True)

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
			clustering_sim=constants[4],
			clustering_cov=constants[5],
			size_threshold = constants[9],
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

		pf.print_message(f"Excluding clusters based on the frequency of the CDSs in the genomes (<={constants[2]})...", "info")
		# Sum the frequency of the CDSs in each cluster to get the frequency per cluster
		frequency_in_genomes: Dict[str, int] = {repid: sum(cds_frequency[seqid] for seqid in members) for repid, members in clusters.items()}
		# Drop clusters with a frequency below the given threshold and add the corresponding CDS IDs to the dropped dict
		infrequent: List[str] = []
		excluded_clusters: List[str] = []
		for repid in clusters:
			if frequency_in_genomes[repid] <= constants[2]:
				infrequent.extend(clusters[repid])
				excluded_clusters.append(repid)
		# Exclude clusters
		for repid in excluded_clusters:
			del clusters[repid]

		# Add seqids excluded based on cluster frequency on genomes
		for seqid in infrequent:
			excluded.setdefault(seqid, []).append("infrequent")

		pf.print_message(f"Excluded {len(excluded_clusters)} clusters (a total of {len(infrequent)} CDSs).")
		pf.print_message(f"{len(clusters)} clusters remain after filtering based on CDS frequency.", "info")

		# Create the files needed to run BLAST
		pf.print_message("Creating input files for BLAST...", "info")
		loci_files, concatenated_files = cof.prepare_cluster_blast_infiles(blast_inputs, clusters, cds_sequences, protein_sequences, blastdb_aliastool_exec)

	if run_mode == 'schema':
		# Check if schemas have loci with the same identifiers
		common_loci = set()
		# Get list of files for the first schema
		schemas_files = [set([ff.file_basename(file, False) for file in ff.get_paths_in_directory_with_suffix(input_schemas[0], ".fasta")])]
		for schema in input_schemas[1:]:
			# Get the list of files for the other schemas and determine if they have loci identifers in common with previous schemas
			current_schema_files = set([ff.file_basename(file, False) for file in ff.get_paths_in_directory_with_suffix(schema, ".fasta")])
			for i, file_set in enumerate(schemas_files):
				common_loci_ids = set.intersection(file_set, current_schema_files)
				if len(common_loci_ids) > 0:
					pf.print_message(f"The input schemas in {input_schemas[i]} and {schema} have the following loci identifiers in common: {', '.join(common_loci_ids)}", "info")
					common_loci.update(common_loci_ids)
			schemas_files.append(current_schema_files)

		if len(common_loci) > 0:
			pf.print_message(f"Found a total of {len(common_loci)} loci identifiers shared between different schemas.", "info")
			pf.print_message("Please make sure that loci in different schemas have different identifiers.", "info")
			sys.exit(0)

		# Calculate the allele call class counts
		input_class_counts = []
		for ad in allelecall_directory:
			class_counts = compute_class_counts(ad)
			input_class_counts.append(class_counts)

		# Merge dictionaries with class counts into a single dictionary with the class counts for all loci
		loci_class_counts = {}
		for cc in input_class_counts:
			loci_class_counts = loci_class_counts | cc

		# Exclude loci based on their frequency on the genomes and the fraction of special classifications
		for locus, counts in loci_class_counts.items():
			if counts[2] <= constants[2]:
				excluded.setdefault(locus, []).append("short")
			if counts[1] >= constants[3]:
				excluded.setdefault(locus, []).append("special")

		frequency_in_genomes = {locus: counts[2] for locus, counts in loci_class_counts.items()}

####################

		# Copy FASTA files from all schemas into the same directory
		temp_schema = ff.join_paths(processed_indata, ["merged_schema"])
		ff.create_directory(temp_schema)
		temp_schema_short = ff.join_paths(temp_schema, ["short"])
		ff.create_directory(temp_schema_short)
		pf.print_message(f"Copying FASTA files from all schemas to {temp_schema_short}...", "info")
		for schema in input_schemas:
			pf.print_message(f'Copying FASTA files from: {schema}', 'info')
			schema_short = ff.join_paths(schema, ["short"])
			# List FASTA files in the "short" directory
			schema_short_files = ff.get_paths_in_directory_with_suffix(schema_short, ".fasta")
			for file in schema_short_files:
				locus_id = ff.file_basename(file).rsplit('_short', 1)[0]
				if locus_id not in excluded:
					shutil.copy(file, temp_schema_short)
				else:
					pf.print_message(f"File {file} was excluded because of the following reasons: {excluded[locus_id]}.", "info")

		pf.print_message("Creating input files for BLAST...", "info")
		loci_files, concatenated_files = cof.prepare_loci_blast_infiles(temp_schema, blast_inputs, constants[6], constants[7], blastdb_aliastool_exec)

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
	self_blastp_inputs = {locus: paths[1] for locus, paths in loci_files.items()}
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
			# Filter results based on criteria

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
			# Filter results based on criteria

			blastn_data = {}
			if blastn_lines:
				# Group all HSPSs for the same query-subject pair
				for line in blastn_lines:
					blastn_data.setdefault((line[0], line[1]), []).append(line[2:])

		# No matches with both BLASTp and BLASTn
		if not blastp_data and not blastn_data:
			pf.print_message(f"No BLASTp and BLASTn matches for {file}.", "info")
			continue

		# Merge BLASTp and BLASTn pairs
		blast_pairs = set.union(set(blastn_data.keys()), set(blastp_data.keys()))

		for qs_pair in blast_pairs:
			query, subject = qs_pair
			blastp_matches = blastp_data.get(qs_pair)
			blastn_matches = blastn_data.get(qs_pair)
			# Determine and add global values
			if qs_pair not in blast_results_dict and (subject, query) not in blast_results_dict:
				# Get query and subject frequency values to the line
				# Use the seqid in the BLAST output files if running in the "cds" mode
				# Determine the locus identifier otherwise
				query_frequency = frequency_in_genomes.get(query, frequency_in_genomes[query.rsplit('_', 1)[0]])
				subject_frequency = frequency_in_genomes.get(subject, frequency_in_genomes[subject.rsplit('_', 1)[0]])

				# Calculate the frequency ratio
				# The frequency of the query or subject should not be 0 because of the frequency threshold initially applied
				frequency_ratio = round(min(query_frequency/subject_frequency, subject_frequency/query_frequency), 3)

				# Get the sequence length for the query and subject
				# Need to get sequence length value from BLASTn matches it there were no BLASTp matches
				if blastp_matches:
					query_len = blastp_matches[0][0]
					subject_len = blastp_matches[0][1]
				else:
					query_len = str(int(int(blastn_matches[0][0])/3))
					subject_len = str(int(int(blastn_matches[0][1])/3))

				global_data = [query_len, subject_len, query_frequency, subject_frequency, frequency_ratio]
				# Create query-subject pair entry
				blast_results_dict.setdefault(qs_pair, [[[],[]], [[],[]], global_data])

			# Get BLASTp data
			if blastp_matches:
				for m in blastp_matches:
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
			if blastn_matches:
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
		# Merge BLASTp aligned intervals
		# Query to subject alignments
		# Query
		query_aligned_intervals: List[int] = sorted([m[1] for m in match_data[0][0]], key= lambda x: x[0])
		# It is possible that there are not BLASTp alignments, only BLASTn
		if len(query_aligned_intervals) > 0:
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
		else:
			blast_results_dict[qs_pair][0].append(None)

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
			blast_results_dict[qs_pair][1].append(None)

		# Merge BLASTn aligned intervals
		# Query to subject alignments
		# Query
		query_aligned_intervals: List[int] = sorted([m[1] for m in match_data[0][1]], key= lambda x: x[0])
		# It is possible that there are not BLASTn alignments
		if len(query_aligned_intervals) > 0:
			query_aligned_merged_intervals = af.merge_intervals(query_aligned_intervals)
			query_aligned_length = sum([(i[1]-i[0]) for i in query_aligned_merged_intervals])
			# Subject
			subject_aligned_intervals: List[int] = sorted([m[2] for m in match_data[0][1]], key= lambda x: x[0])
			subject_aligned_merged_intervals = af.merge_intervals(subject_aligned_intervals)
			subject_aligned_length = sum([(i[1]-i[0]) for i in subject_aligned_merged_intervals])
			# Compute the weighted pident based on the pident of all HSPSs
			# Different HSPSs can overlap and it is not possible to compute the true pident based on BLAST's outmft 6
			query_aligned_identity: float = sum(float(m[3])*((m[1][1]-m[1][0])/(query_length*3)) for m in match_data[0][1])
			subject_aligned_identity: float = sum(float(m[3])*((m[2][1]-m[2][0])/(subject_length*3)) for m in match_data[0][1])
			# Compute the global palign minimum and maximum values
			global_talign_min: float = min(query_aligned_length / query_length,
										subject_aligned_length / subject_length)
			global_talign_max: float = max(query_aligned_length / query_length,
										subject_aligned_length / subject_length)
			# Compute the minimum and maximum global percent identity
			global_pident_min: float = min(query_aligned_identity, subject_aligned_identity)
			global_pident_max: float = max(query_aligned_identity, subject_aligned_identity)
			# Get the maximum BSR value from all query-subject HSPSs
			max_bsr = 0
			blast_results_dict[qs_pair][0].append([round(global_talign_min, 3), round(global_talign_max, 3),
												round(global_pident_min, 3), round(global_pident_max, 3), max_bsr])
		else:
			blast_results_dict[qs_pair][0].append(None)

		# Subject to query alignments
		# Query
		query_aligned_intervals: List[int] = sorted([m[1] for m in match_data[1][1]], key= lambda x: x[0])
		# It is possible that there are no Subject-Query alignments
		if len(query_aligned_intervals) > 0:
			query_aligned_merged_intervals = af.merge_intervals(query_aligned_intervals)
			query_aligned_length = sum([(i[1]-i[0]) for i in query_aligned_merged_intervals])
			# Subject
			subject_aligned_intervals: List[int] = sorted([m[2] for m in match_data[1][1]], key= lambda x: x[0])
			subject_aligned_merged_intervals = af.merge_intervals(subject_aligned_intervals)
			subject_aligned_length = sum([(i[1]-i[0]) for i in subject_aligned_merged_intervals])
			# Compute the weighted pident based on the pident of all HSPSs
			# Different HSPSs can overlap and it is not possible to compute the true pident based on BLAST's outmft 6
			query_aligned_identity: float = sum(float(m[3])*((m[1][1]-m[1][0])/(query_length*3)) for m in match_data[1][1])
			subject_aligned_identity: float = sum(float(m[3])*((m[2][1]-m[2][0])/(subject_length*3)) for m in match_data[1][1])
			# Compute the global palign minimum and maximum values
			global_talign_min: float = min(query_aligned_length / query_length,
										subject_aligned_length / subject_length)
			global_talign_max: float = max(query_aligned_length / query_length,
										subject_aligned_length / subject_length)
			# Compute the minimum and maximum global percent identity
			global_pident_min: float = min(query_aligned_identity, subject_aligned_identity)
			global_pident_max: float = max(query_aligned_identity, subject_aligned_identity)
			# Get the maximum BSR value from all query-subject HSPSs
			max_bsr: float = 0
			blast_results_dict[qs_pair][1].append([round(global_talign_min, 3), round(global_talign_max, 3),
												round(global_pident_min, 3), round(global_pident_max, 3), max_bsr])
		else:
			blast_results_dict[qs_pair][1].append(None)

#########################

	# Assign classes to all matches
	for qs_pair, match_data in blast_results_dict.items():
		query_freq = match_data[2][2]
		subject_freq = match_data[2][3]
		frequency_ratio = match_data[2][4]

		# Get data for Query-Subject matches
		qs_data = match_data[0][2]
		# Get data for Subject-Query matches
		sq_data = match_data[1][2]

		# Classify Query-Subject matches
		# Classify based on global_talign_min and bsr
		if qs_data:
			if qs_data[0] >= constants[0]:
				if qs_data[4] >= constants[8]:
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
		# No match with BLASTp, only with BLASTn
		else:
			assigned_class = '6'

		blast_results_dict[qs_pair][0].append(assigned_class)

		# Classify Subject-Query matches
		# Classify based on global_talign_min and bsr
		if sq_data:
			if sq_data[0] >= constants[0]:
				if sq_data[4] >= constants[8]:
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
		else:
			assigned_class = '6'

		blast_results_dict[qs_pair][1].append(assigned_class)

	# Write output file with information for all matches
	match_data_lines = [ct.MATCHES_HEADER]
	for qs_pair, match_data in blast_results_dict.items():
		query, subject = qs_pair
		global_data = match_data[2]
		query_class = match_data[0][-1]
		subject_class = match_data[1][-1]
		query_blastp_data = match_data[0][-3] if match_data[0][-3] is not None else ["NA"]*5
		subject_blastp_data = match_data[1][-3] if match_data[1][-3] is not None else ["NA"]*5
		query_blastn_data = match_data[0][-2][:-1] if match_data[0][-2] is not None else ["NA"]*4
		subject_blastn_data = match_data[1][-2][:-1] if match_data[1][-2] is not None else ["NA"]*4

		# Merge values into a single list
		qs_pair_data = [query, subject] + list(map(str, global_data)) + [query_class, subject_class] + list(map(str, query_blastp_data)) + list(map(str, subject_blastp_data)) + list(map(str, query_blastn_data)) + list(map(str, subject_blastn_data))
		match_data_lines.append(qs_pair_data)

	match_data_outlines = ["\t".join(line) for line in match_data_lines]
	match_data_outfile = ff.join_paths(output_directory, [ct.MATCHES_FILENAME])
	with open(match_data_outfile, "w") as outfile:
		outfile.write("\n".join(match_data_outlines)+"\n")

#########################

	pf.print_message("Assigning actions based on classifications...", "info")
	assigned_classes = {}
	for qs_pair, match_data in blast_results_dict.items():
		query_class = match_data[0][-1]
		query_frequency = match_data[2][2]
		subject_class = match_data[1][-1]
		subject_frequency = match_data[2][3]
		# Need to group assigned classes for each locus if running in modes "schema" or "schema-vs-schema"
		if run_mode == "cds":
			# When running 
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
		elif top_class in {'1c', '2b', '3b', '4b', '6'}:
			actions[seqid] = ("Choice", top_class)
		# Classes 4c and 5 are Add
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
		# Get the action for the representative CDS of for the locus if running in the "schema" mode
		if run_mode == "cds":
			current_group = [(seqid, *actions[seqid]) for seqid in cc]
		else:
			seqids = set([seqid.rsplit('_', 1)[0] for seqid in cc])
			current_group = [(seqid, *actions[seqid]) for seqid in seqids]
		groups.append(current_group)

	# Add remaining representative seqids as Add and assign class 7
	if run_mode == "cds":
		representative_seqids = clusters.keys()
		# Determine seqids that are not in any component/group
		component_seqids = set()
		for cc in connected_components:
			component_seqids.update(cc)
	# Add loci that had no matches
	else:
		representative_seqids = loci_files.keys()
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

	# Need to create FASTA files with the set of representative alleles listed in the file with recommendations
	###################

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


def main(input_schemas: List[str], output_directory: str, allelecall_directory: List[str],
		annotations: str, alignment_ratio_threshold: float, pident_threshold: float,
		clustering_sim_threshold: float, clustering_cov_threshold:float, genome_presence: int,
		special_classifications: float, absolute_size: int, translation_table: int, bsr: float, size_ratio: float, run_mode: str,
		cpu: int, no_cleanup: bool) -> None:
	"""
	Main function to identify spurious genes in a schema.

	Parameters
	----------
	input_schemas : List[str]
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
							special_classifications,
							clustering_sim_threshold,
							clustering_cov_threshold,
							absolute_size,
							translation_table,
							bsr,
							size_ratio]

	identify_spurious_genes(input_schemas,
							output_directory,
							allelecall_directory,
							annotations,
							constants,
							run_mode,
							cpu,
							no_cleanup)
