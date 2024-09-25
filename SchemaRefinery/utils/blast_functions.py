import subprocess
import os
import sys

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
 
	print(stderr)

	# Exit if it is not possible to create BLAST db
	if len(stderr) > 0:
		sys.exit(f'Error while running BLASTp for {fasta_file}\n'
				 f'{blast_path} returned the following error:\n{stderr}')

	return [stdout, stderr]

def run_blast_with_args_only(blast_args):
    """
    Runs BLAST based on input arguments.

    Parameters
    ----------
    blast_args : list
        Contains list with arguments used in subprocess.

    Returns
    -------
    return : None
        Generates BLAST files at output directory
    """
    
    blast_proc = subprocess.Popen(blast_args,
                        stdout=subprocess.PIPE,
                        stderr=subprocess.PIPE)

    stderr = blast_proc.stderr.readlines()
    if len(stderr) > 0:
        print(stderr)

def run_blast_fastas_multiprocessing(id_, blast_exec, blast_results,
                                                file_dict, all_fasta_file):
    """
    This function, runs blast of representatives of the loci vs consolidation of all of the representatives in single file.

    Parameters
    ----------
    id : str
        id of the locus that will be blasted against all of the representatives sequences.
    blast_exec : str
        Path to the BLAST executable.
    blast_results_all_representatives : str
        Path to the folder were to store blast results.
    representative_file_dict : dict
        Dict that contains the path to file for each sequence(key).
    all_representatives_file : str
        Path to the file of all of the sequences to BLAST against.

    Returns
    -------
    return : list
        locus : str
        blast_results_file : str
        list containing locus id and path to the blast_results_file for that locus.
    """

    blast_results_file = os.path.join(blast_results, f"blast_results_{id_}.tsv")
    
    blast_args = [blast_exec, '-query', file_dict[id_],
                  '-subject',
                  all_fasta_file,
                  '-outfmt',
                  '6 qseqid sseqid qlen slen qstart qend sstart send length score gaps pident',
                  '-out', 
                  blast_results_file]

    run_blast_with_args_only(blast_args)

    return [id_, blast_results_file]

def run_blastdb_multiprocessing(blast_exec, blast_db, fasta_file, id_, blast_output,
              max_hsps=None, threads=1, ids_file=None, blast_task=None,
              max_targets=None):
    """
    Execute BLAST.

    Parameters
    ----------
    blast_path : str
        Path to the BLAST executable.
    blast_db : str
        Path to the BLAST database.
    fasta_file : str
        Path to the Fasta file that contains the sequences
        to align against the database.
    id_ : str
        Identifier of the sequence.
    blast_output : str
        Path to the output file.
    max_hsps : int, optional
        Maximum number of High-Scoring Pairs.
    threads : int, optional
        Number of threads passed to BLAST.
    ids_file : path, optional
        Path to a file with the identifiers of the sequences
        to align against. Used to specify the database sequences
        we want to align against.
    blast_task : str, optional
        BLAST task. Allows to set default parameters for a specific
        type of search.
    max_targets : int, optional
        Maximum number of targets sequences to align against.

    Returns
    -------
    stderr : list
        List with the warnings/errors reported by BLAST.
    """
    blast_results_file = os.path.join(blast_output, f"blast_results_{id_}.tsv")

    blast_args = [blast_exec, '-db', blast_db, '-query', fasta_file,
                '-out', blast_results_file, '-outfmt', 
                '6 qseqid sseqid qlen slen qstart qend sstart send length score gaps pident',
                '-num_threads', str(threads), '-evalue', '0.001']
    if max_hsps is not None:
        blast_args.extend(['-max_hsps', str(max_hsps)])
    if ids_file is not None:
        blast_args.extend(['-seqidlist', ids_file])
    if blast_task is not None:
        blast_args.extend(['-task', blast_task])
    if max_targets is not None:
        blast_args.extend(['-max_target_seqs', str(max_targets)])

    run_blast_with_args_only(blast_args)

    return [id_, blast_results_file]

def run_self_score_multiprocessing(id_, blast_exec, file_path, output):
    blast_results_file = os.path.join(output, f"blast_results_{id_}.tsv")
    
    blast_args = [blast_exec, '-query', file_path,
                  '-subject',
                  file_path,
                  '-outfmt',
                  '6 qseqid sseqid qlen slen qstart qend sstart send length score gaps pident',
                  '-out', 
                  blast_results_file]

    run_blast_with_args_only(blast_args)

    return [id_, blast_results_file]

def compute_bsr(subject_score, query_score):
    """Compute the BLAST Score Ratio for an alignment between two sequences.

    Parameters
    ----------
    subject_score : float
        Alignment raw score computed by BLAST.
    query_score : float
        Raw score computed by BLAST for the
        self-alignment of the query sequence.

    Returns
    -------
    bsr : float
        BLAST Score Ratio for the alignment.
    """
    bsr = subject_score / query_score

    return bsr

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