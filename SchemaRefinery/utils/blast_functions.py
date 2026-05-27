import sys
import subprocess
import concurrent.futures
from typing import Dict, Any, List, Tuple, Union, Optional

try:
	from utils import (file_functions as ff,
					   print_functions as pf,
					   alignments_functions as af)
except ModuleNotFoundError:
	from SchemaRefinery.utils import (file_functions as ff,
								   	  print_functions as pf,
									  alignments_functions as af)


def make_blast_db(makeblastdb_path: str, input_fasta: str, output_path: str, db_type: str) -> Tuple[bytes, Union[bytes, str]]:
	"""
	Create a BLAST database.

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
		Type of the database, nucleotide (nucl) or protein (prot).

	Returns
	-------
	stdout : bytes
		BLAST stdout.
	stderr : bytes or str
		BLAST stderr.
	"""
	# Use '-parse-seqids' to be able to specify sequences to align against
	# Use v5 databases (text file with list of sequence IDs needs to be converted with blastdb_aliastool)
	# Decent performance with all BLAST versions, except v2.11 which runs much slower for unknown reasons
	# BLAST <= 2.11 cannot create v4 databases if sequence IDs are alphanumeric and composed of 4 chars
	# v5 databases accept those IDs but replace '-' with '_', which is an issue when chewie is looking for the original IDs
	makedb_cmd: List[str] = [makeblastdb_path, '-in', input_fasta,
							 '-out', output_path, '-parse_seqids',
							 '-dbtype', db_type, '-blastdb_version', '5']

	makedb_process: subprocess.Popen = subprocess.Popen(makedb_cmd,
													stdout=subprocess.PIPE,
													stderr=subprocess.PIPE)
	stdout: bytes
	stderr: bytes
	stdout, stderr = makedb_process.communicate()

	return stdout, stderr


def run_blastdbcmd(blastdbcmd_path, blast_db, output_file) -> Tuple[bytes, Union[bytes, str]]:
	"""Run blastdbcmd to extract sequences from a BLAST database.

	Parameters
	----------
	blastdbcmd_path : str
		Path to the blastdbcmd executable.
	blast_db : str
		Path to the BLAST database.
	output_file : str
		Path to the output file that will store the sequences.

	Returns
	-------
	stdout : bytes
		BLAST stdout.
	stderr : bytes or str
		BLAST stderr.
	"""
	blastdbcmd_args = [blastdbcmd_path, '-db', blast_db, '-out', output_file, '-entry', 'all']

	blastdbcmd_process = subprocess.Popen(blastdbcmd_args,
										  stdout=subprocess.PIPE,
										  stderr=subprocess.PIPE)

	stdout, stderr = blastdbcmd_process.communicate()

	# Exit if it is not possible to extract sequences from BLAST db
	if len(stderr) > 0:
		sys.exit(f'Cound not extract sequences from {blast_db}.\n')

	return stdout, stderr


def run_blast(blast_args: List[str]) -> None:
	"""
	Runs BLAST based on input arguments.

	Parameters
	----------
	blast_args : list
		Contains list with arguments used in subprocess.

	Returns
	-------
	None
	"""
	blast_proc: subprocess.Popen = subprocess.Popen(blast_args[2],
													stdout=subprocess.PIPE,
													stderr=subprocess.PIPE)

	stderr: List[bytes] = blast_proc.stderr.readlines()
	if len(stderr) > 0:
		pf.print_message(stderr, 'error')

	return [blast_args[0], blast_args[1]]


def run_self_score_multiprocessing(locus_path: Dict[str, str], blast_exec: str, output_folder: str, cpu: int) -> Dict[str, int]:
	"""
	Calculate self-score for each loci using BLASTp.

	Parameters
	----------
	locus_path : dict
		Dictionary with keys as loci identifiers and values as paths to the loci FASTA files.
	blast_exec : str
		Path to the BLASTp executable.
	output_folder : str
		Path to the output folder where results will be stored.
	cpu : int
		Number of CPU cores to use for multiprocessing.

	Returns
	-------
	dict
		Dictionary with loci identifiers as keys and their self-scores as values.
	"""
	blast_results: List[str] = [] # List to store paths to BLAST results
	i: int = 1
	# Create BLAST cmd data for multiprocessing
	additional_args = ['-task', 'blastp-fast', '-qcov_hsp_perc', '100', '-subject_besthit', '-xdrop_ungap', '1', '-xdrop_gap', '1', '-xdrop_gap_final', '1', '-gapopen', '11', '-gapextend', '2']
	blast_input_cmds = [create_blast_cmd(blast_exec, file, locus, output_folder, None, max_hsps=1, additional_args=additional_args) for locus, file in locus_path.items()]

	# Calculate self-score
	with concurrent.futures.ProcessPoolExecutor(max_workers=cpu) as executor:
		for res in executor.map(run_blast, blast_input_cmds):
				blast_results.append(res[1])
				# Print progress
				pf.print_message(f"Running BLAST to calculate self-score for {res[0]}...", "info", end='\r', flush=True)
				i += 1

	# Print to avoid printing next message in the same line as last progress message
	pf.print_message(f"", "info")

	self_scores: Dict[str, int] = {}
	for r in blast_results:
		# Extract self-score from BLAST results
		self_scores.update(af.get_self_scores(r))

	return self_scores


def compute_bsr(subject_score: float, query_score: float) -> float:
	"""
	Compute the BLAST Score Ratio for an alignment between two sequences.

	Parameters
	----------
	subject_score : float
		Alignment raw score computed by BLAST.
	query_score : float
		Raw score computed by BLAST for the
		self-alignment of the query sequence.

	Returns
	-------
	float
		BLAST Score Ratio for the alignment.
	"""
	bsr: float = subject_score / query_score

	return bsr


def determine_blast_task(sequences: List[str], blast_type: str = 'blastp') -> str:
	"""
	Determine the type of BLAST task to execute.

	It is necessary to define the BLAST task if any of the
	sequences to align is shorter than 50 base pairs for
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
	str
		A string that indicates the type of BLAST task to
		execute based on the minimum sequence size.

	Notes
	-----
	More information about the task option at:
		https://www.ncbi.nlm.nih.gov/books/NBK569839/
	"""
	# Get sequence length threshold for BLAST application
	length_threshold: int = 50 if blast_type == 'blastn' else 30
	sequence_lengths: List[int] = [len(p) for p in sequences]
	minimum_length: int = min(sequence_lengths)
	if minimum_length < length_threshold:
		blast_task: str = f'{blast_type}-short'
	else:
		blast_task = blast_type

	return blast_task


def run_blastdb_aliastool(blastdb_aliastool_path: str, seqid_infile: str, seqid_outfile: str) -> Tuple[bytes, Union[bytes, str]]:
	"""
	Convert list of sequence identifiers into binary format.

	Parameters
	----------
	blastdb_aliastool_path : str
		Path to the blastdb_aliastool executable.
	seqid_infile : str
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
	blastdb_aliastool_args: List[str] = [blastdb_aliastool_path, '-seqid_file_in',
										 seqid_infile, '-seqid_file_out', seqid_outfile]

	blastdb_aliastool_process = subprocess.Popen(blastdb_aliastool_args,
												 stdout=subprocess.PIPE,
												 stderr=subprocess.PIPE)

	stdout, stderr = blastdb_aliastool_process.communicate()

	# Exit if it is not possible to create BLAST db
	if len(stderr) > 0:
		sys.exit(f'Could not convert {seqid_infile} to binary format. {blastdb_aliastool_path} returned the following error: {stderr}')

	return stdout, stderr


def run_blast_multiprocessing(input_files: List[str], output_directory: str, blast_db: str, blast_exec: str, cpu: int, max_hsps: int, max_targets: int, additional_args: List[str] = None) -> List[str]:
	"""
	Run BLAST in parallel.

	Parameters
	----------
	input_files : dict
		Dictionary with IDs as keys and list of paths to the fasta file
		and the file with sequence IDs as values.
	output_directory : str
		Path to the folder where to store the BLAST results.
	blast_db : str
		Path to the BLAST database.
	blast_exec : str
		Path to either the BLASTn or BLASTp executable.
	cpu : int
		Number of CPU cores to use for multiprocessing.

	Returns
	-------
	blast_results_files : List[str]
		List containing the paths to the BLAST results files.
	"""
	i: int = 1
	blast_results_files: List[str] = []

	input_ids = list(input_files.keys())
	fasta_files = list([paths[0] for paths in input_files.values()])
	seqids_files = list([paths[1] for paths in input_files.values()])
	# Define the total number of BLASTs for progress message
	total_blasts = len(fasta_files)
	# Determine maximum ID length to print fixed length message (avoids partial ID overlap when previous ID was longer)
	max_id_length = max([len(i) for i in input_ids])

	# Create BLAST inputs
	blast_input_cmds = [create_blast_cmd(blast_exec, fasta_files[i], seqid, output_directory, cpu, blast_db=blast_db, ids_file=seqids_files[i], max_targets=max_targets, max_hsps=max_hsps, additional_args=additional_args) for i, seqid in enumerate(input_ids)]

	with concurrent.futures.ProcessPoolExecutor(max_workers=cpu) as executor:
		for res in executor.map(run_blast, blast_input_cmds):
			# Save the path to the BLAST results file
			blast_results_files.append(res[1])
			pf.print_message(f"Running BLAST for {res[0]} - {i}/{total_blasts:<{max_id_length}}", "info", end='\r', flush=True)
			i += 1

	# Print to avoid printing next message in the same line as last progress message
	pf.print_message(f"", "info")

	return blast_results_files


def create_blast_cmd(blast_exec: str, fasta_file: str, input_id: str, output_directory: str,
					 threads: int, blast_db: str = None, ids_file: str = None, max_targets: int = None,
					 max_hsps: int = None, blast_task: str = None, additional_args: List[str] = None) -> List[str]:
	"""
	Create the command to run BLAST.

	Returns
	-------
	list
		List containing the command to run BLAST.
	"""
	blast_results_file: str = ff.join_paths(output_directory, [f"blast_results_{input_id}.tsv"])

	blast_args: List[str] = [blast_exec,
							 '-query', fasta_file,
							 '-out', blast_results_file,
							 '-outfmt', '6 qseqid sseqid qlen slen qstart qend sstart send length score gaps pident',
							 '-evalue', '0.001']

	if threads:
		blast_args.extend(['-num_threads', str(threads)])

	if blast_db:
		blast_args.extend(['-db', blast_db])
	else:
		blast_args.extend(['-subject', fasta_file])

	if ids_file is not None:
		blast_args.extend(['-negative_seqidlist', ids_file])
	if max_targets is not None:
		blast_args.extend(['-max_target_seqs', str(max_targets)])
	if max_hsps is not None:
		blast_args.extend(['-max_hsps', str(max_hsps)])
	if blast_task is not None:
		blast_args.extend(['-task', blast_task])
	if additional_args is not None:
		blast_args.extend(additional_args)

	return [input_id, blast_results_file, blast_args]
