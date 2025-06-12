import os
from typing import Dict, List, Tuple

try:
	from utils import (
						sequence_functions as sf,
						blast_functions as bf,
						linux_functions as lf,
						file_functions as ff,
						alignments_functions as af,
						Types as tp,
						print_functions as pf,
						logger_functions as logf,
						globals as gb)
except ModuleNotFoundError:
	from SchemaRefinery.utils import (
									sequence_functions as sf,
									blast_functions as bf,
									linux_functions as lf,
									file_functions as ff,
									alignments_functions as af,
									Types as tp,
									print_functions as pf,
									logger_functions as logf,
									globals as gb)


def run_blasts_match_schemas(query_translations_paths: Dict[str, str], blast_db_files: str,
							 blast_folder: str, self_score_dict: Dict[str, float], max_id_length: int,
							 get_blastp_exec: str, bsr: float, cpu: int) -> Dict[str, Dict[str, float]]:
	"""
	Run BLASTp to match schemas and calculate BSR values.

	Parameters
	----------
	query_translations_paths : dict
		Dictionary with keys as query identifiers and values as paths to the query translation files.
	blast_db_files : str
		Path to the BLAST database files.
	blast_folder : str
		Path to the folder where BLAST results will be stored.
	self_score_dict : dict
		Dictionary with query identifiers as keys and their self-scores as values.
	max_id_length : int
		Maximum length of the query identifiers.
	get_blastp_exec : str
		Path to the BLASTp executable.
	bsr : float
		BSR threshold value.
	cpu : int
		Number of CPU cores to use for multiprocessing.

	Returns
	-------
	dict
		Dictionary with loci identifiers as keys and tuples of the best subject identifier and BSR value as values.
	"""
	# Run BLASTp
	blastp_results_folder: str = os.path.join(blast_folder, 'blastp_results')
	ff.create_directory(blastp_results_folder)

	# Initialize dictionaries to store BSR values
	bsr_values: Dict[str, Dict[str, float]] = {}
	best_bsr_values: Dict[str, Dict[str, float]] = {}
	total_blasts: int = len(query_translations_paths)
	blastp_results_files = bf.run_blastp_operations(cpu,
													get_blastp_exec,
													blast_db_files,
													query_translations_paths,
													blastp_results_folder,
													total_blasts,
													max_id_length)

	matched_subjects = set()
	for blast_result_file in blastp_results_files:
		# Get the alignments
		filtered_alignments_dict: tp.BlastDict 
		filtered_alignments_dict, _, _ = af.get_alignments_dict_from_blast_results_simplified(blast_result_file,
																									0,
																									False,
																									True)

		# Since BLAST may find several local alignments, choose the largest one to calculate BSR.
		for query, subjects_dict in filtered_alignments_dict.items():
			# Get the loci name
			query_locus_id: str = query.split('_')[0]

			best_bsr_values.setdefault(query_locus_id, {})
			# Create the dict of the query
			bsr_values.setdefault(query, {})
			for subject_id, results in subjects_dict.items():
				subject_locus_id = '_'.join(subject_id.split('_')[:-1])
				# Highest score (First one)
				subject_score: float = next(iter(results.values()))['score']
				# Calculate BSR value
				computed_score: float = bf.compute_bsr(subject_score, self_score_dict[query])
				# Check if the BSR value is higher than the threshold
				if computed_score >= bsr:
					computed_score = round(computed_score, 3)
					# Due to BLASTp compositional adjustments, some BSR values can be greater than 1
					# Convert those to 1.0
					if computed_score > 1.0:
						computed_score = 1.0
					# Save all of the different matches that this query had and their BSR values
					bsr_values[query].update({subject_id: computed_score})
				else:
					continue
				# Save the best match for each query and subject matches
				subject_locus_id = '_'.join(subject_id.split('_')[:-1])
				if not best_bsr_values[query_locus_id].get(subject_locus_id):
					best_bsr_values[query_locus_id][subject_locus_id] = computed_score
					matched_subjects.add(subject_locus_id)
				elif computed_score > best_bsr_values[query_locus_id][subject_locus_id]:
					best_bsr_values[query_locus_id][subject_locus_id] = computed_score

	# Remove entries for queries without results
	best_bsr_values = {query: match for query, match in best_bsr_values.items() if match}

	return best_bsr_values, matched_subjects


def write_best_matches_to_file(match_data, output_directory, match_type):
	"""Write the best matches to a file.

	Parameters
	----------
	match_data : dict
		Dictionary containing the data for the best matches.
	output_directory : str
		Path to the output directory.
	match_type : str
		Type of matches (e.g., 'hashes_dna', 'hashes_prot', 'reps_vs_reps').

	Returns
	-------
	output_file: str
		Path to the output file.
	"""
	# Define path to output file
	output_file = os.path.join(output_directory, match_type+'_matches.tsv')
	match_lists = []
	for k, v in match_data.items():
		current_matches = [[k, k2, str(v2), match_type] for k2, v2 in v.items()]
		match_lists.extend(current_matches)
	match_lines = ['\t'. join(l) for l in match_lists]
	ff.write_lines(match_lines, output_file)

	return output_file


def translate_files(input_files, output_directory, translation_table):
	"""Translate sequences in FASTA files.

	Parameters
	----------
	input_files : dict
		Dictionary with locus identifiers as keys and paths to the FASTA files as values.
	output_directory : str
		Path to the output directory.
	translation_table : int
		Translation table to use for translation.

	Returns
	-------
	translated_paths : dict
		Dictionary with locus identifiers as keys and paths to the translated FASTA files as values.
	protein_hashes : dict
		Dictionary with protein hashes as keys and locus identifiers as values.
	"""
	# Store paths to FASTA files containing translations
	translated_paths: Dict[str, str] = {}
	# Store protein hashes from all translations
	protein_hashes: Dict[str, str] = {}

	i = 0
	for locus, file_path in input_files.items():
		# Read FASTA file
		records: Dict[str, str] = sf.import_sequences(file_path)
		# Create translation file path
		translated_file_path = os.path.join(output_directory, f"{locus}_translation.fasta")
		# Translate sequences, save untranslated sequences, get the protein hashes and update translation dictionary
		_, prot_hashes, _ = sf.translate_seq_deduplicate(records,
														translated_file_path,
														0,
														translation_table,
														True)
		# Update the query translation and protein dictionaries
		if len(prot_hashes) > 0:
			translated_paths[locus] = translated_file_path
			# Add hashes used for exact matching
			protein_hashes.update(prot_hashes)

		# Increment number of translated files
		i += 1
		pf.print_message(f"Translated {i}/{len(input_files)}", "info", end='\r', flush=True)

	# Add a new line after the progress message
	print()

	return translated_paths, protein_hashes


def match_schemas(first_schema_directory: str, second_schema_directory: str, output_directory: str, bsr: float,
				  translation_table: int, cpu: int, no_cleanup: bool, rep_vs_alleles: bool) -> str:
	"""
	Match schemas between query and subject directories.

	Parameters
	----------
	first_schema_directory : str
		Path to the first schema directory.
	second_schema_directory : str
		Path to the second schema directory.
	output_directory : str
		Path to the output directory.
	bsr : float
		BLAST Score Ratio value.
	translation_table : int
		Genetic code used for translation.
	cpu : int
		Number of CPU cores to use.
	no_cleanup : bool
		If True, temporary files will not be removed.
	rep_vs_alleles: bool
		If True then after the rep vs rep Blast the program will run a second Blast with rep vs alleles.

	Returns
	-------
	None
	"""
	# List to store files with results
	results_files = []
	# A schema files
	a_files: Dict[str, str] = ff.map_basename_to_path(ff.get_paths_in_directory_with_suffix(first_schema_directory, '.fasta'))
	a_files_short: Dict[str, str] = ff.map_basename_to_path(ff.get_paths_in_directory_with_suffix(first_schema_directory+'/short', '.fasta'), remove_suffix='_short')
	# B schema files
	b_files: Dict[str, str] = ff.map_basename_to_path(ff.get_paths_in_directory_with_suffix(second_schema_directory, '.fasta'))
	b_files_short: Dict[str, str] = ff.map_basename_to_path(ff.get_paths_in_directory_with_suffix(second_schema_directory+'/short', '.fasta'), remove_suffix='_short')
	# Choose which schema will be the query (the one with the higher average of alleles per loci)
	total_alleles_a = 0
	for qlocus, fasta_path in a_files.items():
		allele_dict: Dict[str, str] = sf.fetch_fasta_dict(fasta_path, False)
		num_alleles = len(allele_dict)
		total_alleles_a += num_alleles
	# Compute average number of alleles per locus for schema A
	avg_a = total_alleles_a/(len(a_files))

	total_alleles_b = 0
	for qlocus, fasta_path in b_files.items():
		allele_dict: Dict[str, str] = sf.fetch_fasta_dict(fasta_path, False)
		num_alleles = len(allele_dict)
		total_alleles_b += num_alleles
	# Compute average number of alleles per locus for schema B
	avg_b = total_alleles_b/(len(b_files))

	# Create lists with schemas data
	schema_a_data = [first_schema_directory, a_files, a_files_short, total_alleles_a, avg_a]
	schema_b_data = [second_schema_directory, b_files, b_files_short, total_alleles_b, avg_b]

	# Select query and subject
	query_schema_data = schema_a_data if avg_a >= avg_b else schema_b_data
	unmatched_queries = set(query_schema_data[1].keys())
	subject_schema_data = schema_a_data if schema_a_data != query_schema_data else schema_b_data
	unmatched_subjects = set(subject_schema_data[1].keys())

	pf.print_message(f"{query_schema_data[0]} set as Query.", "info")
	pf.print_message(f"{subject_schema_data[0]} set as Subject.", "info")
	pf.print_message(f"Total alleles in Query Schema: {query_schema_data[3]}. Total Loci: {len(query_schema_data[1])}. And an average of {round(query_schema_data[4])} alleles per locus.", "info")
	pf.print_message(f"Total alleles in Subject Schema: {subject_schema_data[3]}. Total Loci: {len(subject_schema_data[1])}. And an average of {round(subject_schema_data[4])} alleles per locus.", "info")

	# Create the output directories
	blast_folder: str = os.path.join(output_directory, 'blast_processing')
	ff.create_directory(blast_folder)
	# Directories for complete schemas
	query_translation_folder: str = os.path.join(blast_folder, 'Query_Translation')
	ff.create_directory(query_translation_folder)
	subject_translation_folder: str = os.path.join(blast_folder, 'Subject_Translation')
	ff.create_directory(subject_translation_folder)
	# Directories for the representatives
	query_translation_rep_folder: str = os.path.join(blast_folder, 'Query_Translation_Rep')
	ff.create_directory(query_translation_rep_folder)
	subject_translation_rep_folder: str = os.path.join(blast_folder, 'Subject_Translation_Rep')
	ff.create_directory(subject_translation_rep_folder)

	# Process query FASTA files for the complete query schema
	pf.print_message('Processing the complete Query FASTA files...', 'info')
	query_hashes: Dict[str, str] = {}
	for qlocus, path in query_schema_data[1].items():
		fasta_dict: Dict[str, str] = sf.import_sequences(path)
		hash_dict: Dict[str, str] = {sf.seq_to_hash(v): k for k, v in fasta_dict.items()}
		# Need to consider that same sequence/hash can be represented more than once
		for k, v in hash_dict.items():
			query_hashes.setdefault(k, []).append(v)

	# Process subject FASTA files for the complete subject schema
	pf.print_message('Processing the complete Subject FASTA files...', 'info')
	subject_hashes: Dict[str, str] = {}
	for slocus, path in subject_schema_data[1].items():
		fasta_dict = sf.import_sequences(path)
		hash_dict: Dict[str, str] = {sf.seq_to_hash(v): k for k, v in fasta_dict.items()}
		# Need to consider that same sequence/hash can be represented more than once
		for k, v in hash_dict.items():
			subject_hashes.setdefault(k, []).append(v)

	# -------------------------------------------------------------------
	# Comparision of the Query and Subject DNA hashes (BSR = 1.0)
	# -------------------------------------------------------------------
	# Prepare best BSR values and query translations
	pf.print_message(f"The query schema has {len(query_hashes)} DNA hashes.", "info")
	pf.print_message(f"The subject schema has {len(subject_hashes)} DNA hashes.", "info")
	pf.print_message("Matching DNA hashes between query and subject schema...", "info")

	# Store query and subject loci that matched and were excluded
	matched_queries = set()
	matched_subjects = set()
	total_matches = 0
	# Find common keys (matching DNA hashes)
	common_keys = set(query_hashes) & set(subject_hashes)
	match_data = {}
	for dna_hash in common_keys:
		# Get query locus ID
		# Possible to have multiple query loci with the same hash
		query_loci = ['_'.join(qlocus.split('_')[:-1]) for qlocus in query_hashes[dna_hash]]
		# Get subject locus ID
		# Possible to have multiple subject loci with the same hash
		subject_loci = ['_'.join(slocus.split('_')[:-1]) for slocus in subject_hashes[dna_hash]]
		# Add matches to the match data
		for query_locus in query_loci:
			match_data.setdefault(query_locus, {})
			# Add each subject locus that matched query locus
			for subject_locus in subject_loci:
				match_data[query_locus].setdefault(subject_locus, 1.0)
				matched_subjects.add(subject_locus)
			matched_queries.add(query_locus)
		# Exclude subject loci that matched
		for subject_locus in subject_loci:
			# Exclude main FASTA file
			subject_schema_data[1].pop(subject_locus, None)
			# Exclude FASTA file with representative alleles
			subject_schema_data[2].pop(subject_locus, None)

	# Keep track of queries and subjects that do not have matches
	unmatched_queries = unmatched_queries - matched_queries
	unmatched_subjects = unmatched_subjects - matched_subjects

	# Write results to the best matches file
	if len(match_data) > 0:
		pf.print_message("Writting results to the output file...", "info")
		best_blast_matches_file = write_best_matches_to_file(match_data, output_directory, 'hashes_dna')
		results_files.append(best_blast_matches_file)
		total_matches += len(match_data)
	else:
		pf.print_message("No matches found based on DNA hashes.", "info")

	# Print out stats
	pf.print_message(f"The DNA hash comparison found {len(common_keys)} matches between alleles.", "info")
	pf.print_message(f"{len(matched_subjects)} subject loci had matches and were excluded.", "info")
	pf.print_message(f"{len(subject_schema_data[1])} subject loci will continue to the next step.", "info")

	# Translate query schema
	pf.print_message("Translating query alleles...", "info")
	_, query_protein_hashes = translate_files(query_schema_data[1], query_translation_folder, translation_table)
	pf.print_message("Translating query representative alleles...", "info")
	query_reps_translated_paths, _ = translate_files(query_schema_data[2], query_translation_rep_folder, translation_table)
	# Translate subject schema
	pf.print_message("Translating subject alleles...", "info")
	subject_translated_paths, subject_protein_hashes = translate_files(subject_schema_data[1], subject_translation_folder, translation_table)
	pf.print_message("Translating subject representative alleles...", "info")
	subject_reps_translated_paths, _ = translate_files(subject_schema_data[2], subject_translation_rep_folder, translation_table)

	# -------------------------------------------------------------------
	# Comparision of the Query and Subject protein hashes (the BSR = 1.0)
	# -------------------------------------------------------------------
	# Prepare best BSR values and query translations
	pf.print_message("Matching protein hashes between query and subject schema...", "info")
	pf.print_message(f"The query schema has {len(query_protein_hashes)} protein hashes.", "info")
	pf.print_message(f"The subject schema has {len(subject_protein_hashes)} protein hashes.", "info")
	best_bsr_values = {}

	# Store query and subject loci that matched and were excluded
	matched_queries = set()
	matched_subjects = set()
	# Find common keys (matching protein hashes)
	common_keys = set(query_protein_hashes) & set(subject_protein_hashes)
	match_data = {}
	for prot_hash in common_keys:
		# Get query locus ID
		# Possible to have multiple query loci with the same hash
		query_loci = ['_'.join(qlocus.split('_')[:-1]) for qlocus in query_protein_hashes[prot_hash]]
		# Get subject locus ID
		# Possible to have multiple subject loci with the same hash
		subject_loci = ['_'.join(slocus.split('_')[:-1]) for slocus in subject_protein_hashes[prot_hash]]
		# Add matches to the match data
		for query_locus in query_loci:
			match_data.setdefault(query_locus, {})
			# Add each subject locus that matched query locus
			for subject_locus in subject_loci:
				match_data[query_locus].setdefault(subject_locus, 1.0)
				matched_subjects.add(subject_locus)
			matched_queries.add(query_locus)
		# Exclude subject loci that matched
		for subject_locus in subject_loci:
			# Exclude main translated FASTA file
			subject_translated_paths.pop(subject_locus, None)
			# Exclude FASTA file with representative alleles
			subject_reps_translated_paths.pop(subject_locus, None)

	# Keep track of queries and subjects that do not have matches
	unmatched_queries = unmatched_queries - matched_queries
	unmatched_subjects = unmatched_subjects - matched_subjects

	# Write results to the best matches file
	if len(match_data) > 0:
		pf.print_message("Writting results to the output file...", "info")
		best_blast_matches_file = write_best_matches_to_file(match_data, output_directory, 'hashes_prot')
		results_files.append(best_blast_matches_file)
		total_matches += len(match_data)
	else:
		pf.print_message("No matches found based on protein hashes.", "info")

	# Print out stats
	pf.print_message(f"The protein hash comparison found {len(common_keys)} matches between alleles.", "info")
	pf.print_message(f"{len(matched_subjects)} subject loci had matches and were excluded.", "info")
	pf.print_message(f"{len(subject_translated_paths)} subject loci will continue to the next step.", "info")

	# -------------------------------------------------------------------
	# Blast with rep vs rep
	# -------------------------------------------------------------------
	# Get Path to the blastp executable
	get_blastp_exec: str = lf.get_tool_path('blastp')

	# Get the maximum length of the IDs for better prints
	max_id_length: int = len(max(query_reps_translated_paths.keys(), key=len))

	# Calculate self-scores for each query
	pf.print_message("Calculating self-scores for each locus...", 'info')
	self_score_dict: Dict[str, float] = bf.calculate_self_score(query_reps_translated_paths,
															  get_blastp_exec,
															  blast_folder,
															  max_id_length,
															  cpu)

	# Create BLAST database
	pf.print_message("Creating BLASTp database for subject representatives...", "info")
	# Concatenate FASTA files with subject representatives
	concat_file = os.path.join(subject_translation_rep_folder, 'subject_reps_vs_reps_concat.fasta')
	concat_file = ff.concatenate_files(subject_reps_translated_paths.values(), concat_file, header=None)

	blastdb_path: str = os.path.join(blast_folder, 'subject_reps_vs_reps_blastdb')
	ff.create_directory(blastdb_path)
	blast_db_files: str = os.path.join(blastdb_path, 'subject_reps_vs_reps_proteins_db')
	makeblastdb_exec: str = lf.get_tool_path('makeblastdb')
	bf.make_blast_db(makeblastdb_exec, concat_file, blast_db_files, 'prot')

	# Run BLAST rep_vs_rep
	pf.print_message("Running reps vs reps BLASTp...", "info")
	best_bsr_values, matched_subjects = run_blasts_match_schemas(query_reps_translated_paths,
																			blast_db_files,
																			blast_folder,
																			self_score_dict,
																			max_id_length,
																			get_blastp_exec,
																			bsr,
																			cpu)

	# Write the best BLAST matches to a file
	if len(best_bsr_values) > 0:
		pf.print_message("Writting results to the output file...", "info")
		best_blast_matches_file = write_best_matches_to_file(best_bsr_values, output_directory, 'reps_vs_reps')
		results_files.append(best_blast_matches_file)
		total_matches += len(best_bsr_values)
	else:
		pf.print_message("No matches found based on reps vs reps BLASTp.", "info")

	# Remove matched id from the full subject schema file
	for subject_locus in matched_subjects:
		subject_translated_paths.pop(subject_locus, None)
		subject_reps_translated_paths.pop(subject_locus, None)

	matched_queries = set(best_bsr_values.keys())
	# Keep track of queries and subjects that do not have matches
	unmatched_queries = unmatched_queries - matched_queries
	unmatched_subjects = unmatched_subjects - matched_subjects

	pf.print_message(f"The BLASTp rep vs rep comparison found {len(best_bsr_values)} matches.", "info")
	pf.print_message(f"{len(matched_subjects)} subject loci had matches and were excluded.", "info")
	if rep_vs_alleles:
		pf.print_message(f"{len(subject_reps_translated_paths)} subject loci will continue to the next step.", "info")

	# -------------------------------------------------------------------
	# Blast with rep vs alleles
	# -------------------------------------------------------------------
	# Create BLAST database
	if rep_vs_alleles and len(subject_translated_paths) > 0:
		# Concatenate FASTA files with subject representatives
		concat_file = os.path.join(subject_translation_rep_folder, 'subject_reps_vs_alleles_concat.fasta')
		concat_file = ff.concatenate_files(subject_reps_translated_paths.values(), concat_file, header=None)
		pf.print_message("Creating BLAST database with complete subject schema...", "info")
		blastdb_path: str = os.path.join(blast_folder, 'subject_reps_vs_alleles_blastdb')
		ff.create_directory(blastdb_path)
		blast_db_files: str = os.path.join(blastdb_path, 'subject_reps_vs_alleles_proteins_db')
		makeblastdb_exec: str = lf.get_tool_path('makeblastdb')
		bf.make_blast_db(makeblastdb_exec, concat_file, blast_db_files, 'prot')

		# Run BLAST rep query vs alleles subject
		pf.print_message("Running reps vs alleles BLASTp...", "info")
		best_bsr_values, matched_subjects = run_blasts_match_schemas(query_reps_translated_paths,
																				blast_db_files,
																				blast_folder,
																				self_score_dict,
																				max_id_length,
																				get_blastp_exec,
																				bsr,
																				cpu)

		# Write the best BLAST matches to a file
		if len(best_bsr_values) > 0:
			pf.print_message("Writting results to the output file...", "info")
			best_blast_matches_file = write_best_matches_to_file(best_bsr_values, output_directory, 'reps_vs_alleles')
			results_files.append(best_blast_matches_file)
			total_matches += len(best_bsr_values)

		else:
			pf.print_message("No matches found based on reps vs alleles BLASTp.", "info")

		matched_queries = set(best_bsr_values.keys())
		# Keep track of queries and subjects that do not have matches
		unmatched_queries = unmatched_queries - matched_queries
		unmatched_subjects = unmatched_subjects - matched_subjects

		pf.print_message(f"The BLASTp reps vs alleles comparison found {len(best_bsr_values)} matches.", "info")
		pf.print_message(f"{len(matched_subjects)} subject loci had matches and were excluded.", "info")

	pf.print_message(f"Writing file with list of unmatched queries and subjects...", "info")
	# Create file with list of queries and subjects that had no matches
	unmatched_queries_lines = [f'{query}\tNot matched\tNA\treps_vs_reps'
							   if not rep_vs_alleles
							   else f'{query}\tNot matched\tNA\treps_vs_alleles'
							   for query in unmatched_queries]
	unmatched_subjects_lines = [f'Not matched\t{subject}\tNA\treps_vs_reps'
							   if not rep_vs_alleles
							   else f'Not matched\t{subject}\tNA\treps_vs_alleles'
							   for subject in unmatched_subjects]
	unmatched_lines = unmatched_queries_lines + unmatched_subjects_lines
	if len(unmatched_lines) > 0: 
		unmatched_file = os.path.join(output_directory, 'unmatched.tsv')
		ff.write_lines(unmatched_lines, unmatched_file)
		results_files.append(unmatched_file)

	# Concatenate results files
	output_d= os.path.abspath(output_directory)
	final_results_file = os.path.join(output_d, 'Match_Schemas_Results.tsv')
	ff.concatenate_files(results_files, final_results_file, header='Query\tSubject\tBSR\tProcess\n')

	# Print final statistics
	pf.print_message('')
	pf.print_message(f'A total of {total_matches} loci were matched from each schema.', 'info')
	pf.print_message(f'\t {len(unmatched_queries)} Query loci were not matched.', 'info')
	pf.print_message(f'\t {len(unmatched_subjects)} Subject loci were not matched.', 'info')
	pf.print_message('')

	# Clean up temporary files
	if not no_cleanup:
		pf.print_message("Cleaning up temporary files...", "info")
		# Remove temporary files
		ff.cleanup(output_d, [final_results_file, logf.get_log_file_path(gb.LOGGER)])

	return final_results_file
