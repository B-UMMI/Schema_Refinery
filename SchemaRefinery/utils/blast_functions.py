import subprocess
import os
import sys
from itertools import repeat
import concurrent.futures
from typing import Dict, Any, List, Tuple, Union, Optional


from SchemaRefinery.utils import (file_functions as ff,
                                    blast_functions as bf,
                                    alignments_functions as af,
                                    print_functions as pf)

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

def run_blast(blast_path: str, blast_db: str, fasta_file: str, blast_output: str,
              max_hsps: int = 1, threads: int = 1, ids_file: Optional[str] = None,
              blast_task: Optional[str] = None, max_targets: Optional[int] = None,
              composition_stats: Optional[int] = None) -> Tuple[bytes, Union[bytes, str]]:
    """
    Execute BLAST to align sequences against a BLAST database.

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
    max_hsps : int, optional
        Maximum number of High Scoring Pairs per pair of aligned
        sequences.
    threads : int, optional
        Number of threads/cores used to run BLAST.
    ids_file : str, optional
        Path to a file with sequence identifiers, one per line.
        Sequences will only be aligned to the sequences in the
        BLAST database that match any of the identifiers in this
        file.
    blast_task : str, optional
        Type of BLAST task.
    max_targets : int, optional
        Maximum number of target/subject sequences to align
        against.
    composition_stats : int, optional
        Specify the composition-based statistics method used
        by BLAST.

    Returns
    -------
    stdout : bytes
        BLAST stdout.
    stderr : bytes or str
        BLAST stderr.
    """
    # Do not retrieve hits with high probability of occurring by chance
    blast_args: List[str] = [blast_path, '-db', blast_db, '-query', fasta_file,
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

    blast_process: subprocess.Popen = subprocess.Popen(blast_args,
                                                        stdout=subprocess.PIPE,
                                                        stderr=subprocess.PIPE)

    stdout: bytes
    stderr: bytes
    stdout, stderr = blast_process.communicate()

    pf.print_message(stderr, 'error')

    # Exit if it is not possible to create BLAST db
    if len(stderr) > 0:
        pf.print_message(f'Error while running BLASTp for {fasta_file}, {blast_path} returned the following error: {stderr.decode("utf-8")}', 'error')
        sys.exit()

    return stdout, stderr

def run_blast_with_args_only(blast_args: List[str]) -> None:
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
    
    blast_proc: subprocess.Popen = subprocess.Popen(blast_args,
                                                stdout=subprocess.PIPE,
                                                stderr=subprocess.PIPE)

    stderr: List[bytes] = blast_proc.stderr.readlines()
    if len(stderr) > 0:
        pf.print_message(stderr, 'error')


def run_blast_fastas_multiprocessing(id_: str, blast_exec: str, blast_results: str,
                                     file_dict: Dict[str, str], all_fasta_file: str) -> List[str]:
    """
    Runs BLAST of representatives of the loci vs consolidation of all of the representatives in a single file.

    Parameters
    ----------
    id_ : str
        ID of the locus that will be blasted against all of the representatives sequences.
    blast_exec : str
        Path to the BLAST executable.
    blast_results : str
        Path to the folder where to store BLAST results.
    file_dict : dict
        Dictionary that contains the path to file for each sequence (key).
    all_fasta_file : str
        Path to the file of all of the sequences to BLAST against.

    Returns
    -------
    list
        List containing locus ID and path to the BLAST results file for that locus.
    """

    blast_results_file: str = os.path.join(blast_results, f"blast_results_{id_}.tsv")
    
    blast_args: List[str] = [blast_exec, '-query', file_dict[id_],
                             '-subject', all_fasta_file,
                             '-outfmt', '6 qseqid sseqid qlen slen qstart qend sstart send length score gaps pident',
                             '-out', blast_results_file]

    run_blast_with_args_only(blast_args)

    return [id_, blast_results_file]


def run_blastdb_multiprocessing(blast_exec: str, blast_db: str, fasta_file: str, id_: str, blast_output: str,
                                max_hsps: Optional[int] = None, threads: int = 1, ids_file: Optional[str] = None,
                                blast_task: Optional[str] = None, max_targets: Optional[int] = None) -> List[str]:
    """
    Execute BLAST.

    Parameters
    ----------
    blast_exec : str
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
    ids_file : str, optional
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
    list
        List containing the sequence identifier and the path to the BLAST results file.
    """
    blast_results_file: str = os.path.join(blast_output, f"blast_results_{id_}.tsv")

    blast_args: List[str] = [blast_exec, '-db', blast_db, '-query', fasta_file,
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


def run_self_score_multiprocessing(id_: str, blast_exec: str, file_path: str, output: str) -> List[str]:
    """
    Execute BLAST to calculate self-score.

    Parameters
    ----------
    id_ : str
        Identifier of the sequence.
    blast_exec : str
        Path to the BLAST executable.
    file_path : str
        Path to the Fasta file that contains the sequence.
    output : str
        Path to the output directory.

    Returns
    -------
    list
        List containing the sequence identifier and the path to the BLAST results file.
    """
    blast_results_file: str = os.path.join(output, f"blast_results_{id_}.tsv")
    
    blast_args: List[str] = [blast_exec, '-query', file_path,
                             '-subject', file_path,
                             '-outfmt', '6 qseqid sseqid qlen slen qstart qend sstart send length score gaps pident',
                             '-out', blast_results_file]

    run_blast_with_args_only(blast_args)

    return [id_, blast_results_file]


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


def calculate_self_score(paths_dict: Dict[str, str], blast_exec: str, output_folder: str, max_id_length: int, cpu: int) -> Dict[str, int]:
    """
    Calculate self-score for each loci using BLASTp.

    Parameters
    ----------
    paths_dict : dict
        Dictionary with keys as loci identifiers and values as paths to the loci files.
    blast_exec : str
        Path to the BLASTp executable.
    output_folder : str
        Path to the output folder where results will be stored.
    max_id_length : int
        Maximum length of the loci identifiers.
    cpu : int
        Number of CPU cores to use for multiprocessing.

    Returns
    -------
    dict
        Dictionary with loci identifiers as keys and their self-scores as values.
    """
    
    # Self-score folder
    self_score_folder: str = os.path.join(output_folder, 'self_score_folder')
    ff.create_directory(self_score_folder)

    self_score_dict: Dict[str, Any] = {}
    self_score_results_files: List[str] = [] # List to store paths to BLAST results
    i: int = 1
    # Calculate self-score
    pf.print_message("Calculating self-score for each loci:", 'info')
    with concurrent.futures.ProcessPoolExecutor(max_workers=cpu) as executor:
        for res in executor.map(bf.run_self_score_multiprocessing,
                                paths_dict.keys(),
                                repeat(blast_exec),
                                paths_dict.values(),
                                repeat(self_score_folder)):
            self_score_results_files.append(res[1])
                            
            # Print progress
            pf.print_message(f"Running BLASTp to calculate self-score for {res[0]: <{max_id_length}}", "info", end='\r', flush=True)
            i += 1    

    for blast_results_file in self_score_results_files:
        # Extract self-score from BLAST results
        _, self_score, _, _ = af.get_alignments_dict_from_blast_results(blast_results_file, 0, False, True, True, True, True)

        # Save self-score
        self_score_dict.update(self_score)
    
    pf.print_message("", None)

    return self_score_dict

def run_blastn_operations(cpu: int, get_blastn_exec: str, blast_db: str, rep_paths_nuc: Dict[str, str],
                          all_alleles: Dict[str, List[str]], blastn_results_folder: str,
                          total_reps: int, max_id_length: int) -> List[str]:
    """
    Run BLASTn in parallel for all the cluster representatives.

    Parameters
    ----------
    cpu : int
        Number of CPU cores to use for multiprocessing.
    get_blastn_exec : str
        Path to the BLASTn executable.
    blast_db : str
        Path to the BLAST database files.
    rep_paths_nuc : Dict[str, str]
        Dictionary with the paths to the nucleotide sequences.
    all_alleles : Dict[str, List[str]]
        Dictionary with the alleles for each locus.
    blastn_results_folder : str
        Path to the folder where to store the BLASTn results.
    total_reps : int
        Total number of BLASTn runs to perform.
    max_id_length : int
        Maximum length of the cluster identifiers.
    
    Returns
    -------
    blastn_results_files : List[str]
        List containing the paths to the BLASTn results files.
    """
    blastn_results_files: List[str] = []  # List to store the results files
    i: int = 1
    rep_paths_nuc_list = list(rep_paths_nuc.values())  # Convert dict_values to list
    with concurrent.futures.ProcessPoolExecutor(max_workers=cpu) as executor:
        for res in executor.map(bf.run_blastdb_multiprocessing,
                                repeat(get_blastn_exec),
                                repeat(blast_db),
                                rep_paths_nuc_list,
                                all_alleles,
                                repeat(blastn_results_folder)):
            # Append the results file to the list
            blastn_results_files.append(res[1])
            pf.print_message(f"Running BLASTn for cluster representatives: {res[0]} - {i}/{total_reps: <{max_id_length}}", "info", end='\r', flush=True)
            i += 1
    return blastn_results_files

def run_blastp_operations_based_on_blastn(cpu: int, blastp_runs_to_do, get_blastp_exec: str,
                                          blastp_results_folder: str, rep_paths_prot, rep_matches_prot,
                                          total_blasts: int, max_id_length: int) -> List[str]:
    """
    Run BLASTp in parallel based on the BLASTn results.

    Parameters
    ----------
    cpu : int
        Number of CPU cores to use for multiprocessing.
    blastp_runs_to_do : list
        List with the BLASTp runs to perform.
    get_blastp_exec : str
        Path to the BLASTp executable.
    blastp_results_folder : str
        Path to the folder where to store the BLASTp results.
    rep_paths_prot : dict
        Dictionary with the paths to the protein sequences.
    rep_matches_prot : dict
        Dictionary with the protein sequences that match the nucleotide sequences.
    total_blasts : int
        Total number of BLASTp runs to perform.
    max_id_length : int
        Maximum length of the cluster identifiers.
    
    Returns
    -------
    blastp_results_files : List[str]
        List containing the paths to the BLASTp results files.
    """
    blastp_results_files: List[str] = []  # To store the results files
    i = 1
    rep_matches_prot_list = list(rep_matches_prot.values())  # Convert dict_values to list
    with concurrent.futures.ProcessPoolExecutor(max_workers=cpu) as executor:
        for res in executor.map(bf.run_blast_fastas_multiprocessing,
                                blastp_runs_to_do, 
                                repeat(get_blastp_exec),
                                repeat(blastp_results_folder),
                                repeat(rep_paths_prot),
                                rep_matches_prot_list):
            # Append the results file to the list
            blastp_results_files.append(res[1])
            pf.print_message(f"Running BLASTp for cluster representatives matches: {res[0]} - {i}/{total_blasts: <{max_id_length}}", "info", end='\r', flush=True)
            i += 1
    return blastp_results_files

def run_blastp_operations(cpu: int, get_blastp_exec: str, blast_db_files: str, translations_paths,
                          blastp_results_folder: str, total_blasts: int, max_id_length: int) -> List[str]:
    """
    Run BLASTp in parallel for all the cluster representatives.

    Parameters
    ----------
    cpu : int
        Number of CPU cores to use for multiprocessing.
    get_blastp_exec : str
        Path to the BLASTp executable.
    blast_db_files : str
        Path to the BLAST database files.
    translations_paths : dict
        Dictionary with the paths to the translated sequences.
    blastp_results_folder : str
        Path to the folder where to store the BLASTp results.
    total_blasts : int
        Total number of BLASTp runs to perform.
    max_id_length : int
        Maximum length of the cluster identifiers.
    
    Returns
    -------
    blastp_results_files : List[str]
        List containing the paths to the BLASTp results files.
    """
    blastp_results_files: List[str] = []  # List to store the paths to the BLASTp results files
    i: int = 1
    translations_paths_values = list(translations_paths.values())  # Convert dict_values to list
    translations_paths_keys = list(translations_paths.keys())  # Convert dict_keys to list
    
    with concurrent.futures.ProcessPoolExecutor(max_workers=cpu) as executor:
        for res in executor.map(bf.run_blastdb_multiprocessing,
                                repeat(get_blastp_exec),
                                repeat(blast_db_files),
                                translations_paths_values,
                                translations_paths_keys,
                                repeat(blastp_results_folder)):
            # Save the path to the BLASTp results file
            blastp_results_files.append(res[1])
            pf.print_message(f"Running BLASTp for cluster representatives matches: {res[0]} - {i}/{total_blasts:<{max_id_length}}", "info", end='\r', flush=True)
            i += 1

    return blastp_results_files