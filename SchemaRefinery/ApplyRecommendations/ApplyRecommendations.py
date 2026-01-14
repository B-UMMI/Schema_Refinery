import os
import sys
import csv
import shutil
from typing import Dict, List

try:
	from AdaptLoci import AdaptLoci
	from utils import (sequence_functions as sf,
					   file_functions as ff,
					   logger_functions as logf,
					   print_functions as pf,
					   globals as gb)
except ModuleNotFoundError:
	from SchemaRefinery.AdaptLoci import AdaptLoci
	from SchemaRefinery.utils import (sequence_functions as sf,
									  file_functions as ff,
									  logger_functions as logf,
									  print_functions as pf,
									  globals as gb)


def process_recommendations(recommendations_file: str, 
						    fastas_folder: str,
						    output_directory: str,
						    training_file: str,
						    cpu: int,
						    bsr: float,
						    translation_table: int,
						    no_cleanup:bool) -> None:
	"""
	Creates a schema structure based on the recommendations provided in the recommendations file.

	Parameters
	----------
	recommendations_file : str
		Path to the TSV file containing the recommendations.
	fastas_folder : str
		Path to the folder containing the loci FASTA files.
	output_directory: str
		Path to the directory where the final schema will be stored.
	training_file: str
		Path to the Prodigal training file that will be included in the directory of the adapted schema.
	cpu : int
		Number of CPU cores that will be used to run the process.
	bsr : float
		The BLAST Score Ratio value that will be passed to chewBBACA's PrepExternalSchema module.
	translation_table : int
		Genetic code used to translate alleles.
	no_cleanup : bool
		Flag to indicate whether to clean up temporary files.

	Returns
	-------
	None
		The function writes the output files to the specified directory.
	"""
	# Determine absolute path to output directory (just in case)
	output_directory = os.path.abspath(output_directory)

	# Get all FASTA paths in the FASTA folder
	fasta_files: Dict[str, str] = {
		os.path.basename(fasta_file).split('.')[0]: os.path.join(fastas_folder, fasta_file)
		for fasta_file in os.listdir(fastas_folder) if fasta_file.endswith('.fasta')
	}
	pf.print_message(f'The input schema ({fastas_folder}) has {len(fasta_files)} loci.', 'info')

	temp_fasta_folder = os.path.join(output_directory, 'temp_fastas')
	ff.create_directory(temp_fasta_folder)
	# Read the recommendations file
	with open(recommendations_file, 'r') as f:
		# Skip the header
		recommendations_lines = list(csv.reader(f, delimiter='\t'))[1:]

	# Check if recommendations include Choice actions
	if any([line[1] == 'Choice' for line in recommendations_lines]):
		pf.print_message('The input recommendation file still has loci labeled "Choice".', 'warning')
		pf.print_message('Please change these into Add, Join or Drop.', 'warning')
		sys.exit()

	# Identify each group by splitting on lines that contain only '#'
	group_id: int = 1
	# Initialize the action list dictionary with the first group
	action_list = {group_id: {'Join': [], 'Add': [], 'Drop': []}}
	for line in recommendations_lines:
		if '#' in line:
			group_id += 1
			action_list.setdefault(group_id, {'Join': [], 'Add': [], 'Drop': []})
			continue
		locus, recommendation, class_type = line
		if recommendation not in ['Join', 'Add', 'Drop']:
			pf.print_message(f'Recommendation {recommendation} for locus {locus} (Group {group_id}) is not valid. Ignoring this line.', 'warning')
			continue
		# Add recommendation to the action list if it is valid
		action_list[group_id][recommendation].append(locus)

	# Process each group to generate new FASTA files based on recommendations
	total_drop = 0
	drop_ids: List[str] = []
	total_add = 0
	add_ids: List[str] = []
	total_join = 0
	join_ids: List[str] = []
	for gid, recommendations in action_list.items():
		pf.print_message(f'Processing group {gid}...', 'info')
		pf.print_message(f'{len(recommendations["Add"])} loci to add, {len(recommendations["Drop"])} loci to drop, and {len(recommendations["Join"])} loci to join', 'info')
		for recommendation_type in recommendations:
			current_loci = recommendations[recommendation_type]
			if len(current_loci) > 0:
				pf.print_message(f'Processing {recommendation_type} for loci: {current_loci} in group {gid}', 'info', stdout=False)
				if recommendation_type == 'Drop':
					total_drop += len(current_loci)
					drop_ids.extend(current_loci)
					pf.print_message(f"The following loci will not be included in the schema: {', '.join(current_loci)}.", "info", stdout=False)
				elif recommendation_type == 'Add':
					for locus in current_loci:
						# Define path to FASTA file in the input folder
						fasta_file = fasta_files[locus]
						# Define path to the new FASTA file
						destination_file = os.path.join(temp_fasta_folder, f'{locus}.fasta')
						# Copy the FASTA file to the temporary folder
						shutil.copy(fasta_file, destination_file)
						pf.print_message(f'File {fasta_file} copied to {destination_file}', "info", stdout=False)
						total_add += 1
					add_ids.extend(current_loci)
				elif recommendation_type == 'Join':
					total_join += len(current_loci)
					merged_locusID: str = current_loci[0]
					pf.print_message(f'The following loci FASTA files will be merged under the locus name {current_loci[0]}', "info", stdout=False)
					destination_file: str = os.path.join(temp_fasta_folder, f'{current_loci[0]}.fasta')
					# Initialize the IDs to assign to the alleles in the merged FASTA file
					allele_id: int = 1
					# Initiallize the dictionary to store the merged records
					merged_records: Dict[str, tuple] = {}
					# Count the total number of alleles processed
					total_alleles: int = 0
					for locus in current_loci:
						# The locus must match a FASTA file in the input directory
						if locus in fasta_files:
							# Get the path to the FASTA file
							fasta_file: str = fasta_files[locus]
							# Read the FASTA file and fetch the sequence records
							records: Dict[str, str] = sf.fetch_fasta_dict(fasta_file, False)
							# Add the records to the merged locus
							for seqid, sequence in records.items():
								total_alleles += 1
								# Hash the sequence to get a unique identifier
								sequence_hash: str = sf.hash_sequence(sequence)
								if sequence_hash not in merged_records:
									merged_records[sequence_hash] = (sequence, allele_id)
									# Increment allele_id for the next unique sequence
									allele_id += 1
								else:
									pf.print_message(f'Sequence {seqid} matches a previous sequence already added.', "info", stdout=False)
									continue
						else:
							pf.print_message(f'File {locus} not found in the FASTA folder', "warning")

					# Write the merged records to the new FASTA file
					with open(destination_file, 'w') as outfile:
						outrecords = [f'>{merged_locusID}_{allele_id}\n{sequence}' for sequence_hash, (sequence, allele_id) in merged_records.items()]
						outfile.write('\n'.join(outrecords)+'\n')
					join_ids.extend(current_loci)
			else:
				pf.print_message(f'No loci to process for recommendation {recommendation_type} in group {gid}', 'info', stdout=False)

	# Create schema structure
	pf.print_message("Creating schema structure...", "info")
	# Schema path
	schema_directory = os.path.join(output_directory, 'schema')
	AdaptLoci.adapt_loci(temp_fasta_folder, schema_directory, training_file, cpu, bsr, translation_table)

	# Print total for each action type
	pf.print_message('Summary of actions taken:', 'info')
	pf.print_message(f'\tTotal loci joined: {total_join}', 'info')
	pf.print_message(f'\tTotal loci added: {total_add}', 'info')
	pf.print_message(f'\tTotal loci dropped: {total_drop}', 'info')

	# Determine if there are any loci in the original schema that were not included in the recommendations
	all_ids: List[str] = drop_ids + add_ids + join_ids
	extra_schema_ids: List[str] = [locus_id for locus_id in fasta_files if locus_id not in all_ids]
	if len(extra_schema_ids) != 0:
		pf.print_message(f'{len(extra_schema_ids)} loci in the original schema were not included in the recommendations and were not added to the new schema.', 'info')
		pf.print_message(f'These loci were:\n' + "\n".join(extra_schema_ids), 'info')

	# Count the final number of loci in the schema
	total_loci = [file for file in os.listdir(schema_directory+"/adapted_schema") if file.endswith('.fasta')]
	pf.print_message(f'The final schema has a total of {len(total_loci)} loci.', 'info')

	if not no_cleanup:
		pf.print_message("Cleaning up temporary files...", "info")
		ff.cleanup(output_directory, [schema_directory, logf.get_log_file_path(gb.LOGGER)])

	pf.print_message("ApplyRecommendations process completed.", "info")
