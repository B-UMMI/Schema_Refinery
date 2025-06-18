#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Purpose
-------
This module selects representative alleles for a set of loci.

Code documentation
------------------
"""
import subprocess
import os
import sys
import pathlib
import shutil
from Bio import SeqIO
from typing import Dict, List, Tuple

try:
    from utils import (print_functions as pf,
                        constants as ct,
                        file_functions as ff,
                        blast_functions as bf)
except ModuleNotFoundError:
    from SchemaRefinery.utils import (print_functions as pf,
                                        constants as ct,
                                        file_functions as ff,
                                        blast_functions as bf)

def get_file_prefixes(path_list) -> Dict[str, str]:
	"""
    Determine the file prefix for each file in a list of file paths.

	Parameters
	----------
	path_list : list
		List with file paths.

	Returns
	-------
	prefixes :  dict
		Dictionary with file prefixes as keys and file basenames
		as values.
	"""
	basenames = [pathlib.Path(file).stem for file in path_list]
	prefixes = {}
	for name in basenames:
		prefix = name.rsplit('-protein', 1)[0]
		prefixes.setdefault(prefix, []).append(name)
		#pf.print_message(f'Prefixes: {prefix}', 'info')

	return prefixes


def check_prefix_pdb(input_list, output_directory, makeblastdb_path, blastdbcmd_path) -> None:
	"""
    Check if the BLAST database includes the expected sequence IDs.

	Parameters
	----------
	input_list : str
		Path to file that contains the list of paths to input files.
	output_directory : str
		Path to the directory where dummy data will be created.
	makeblastdb_path : str
		Path to the makeblastdb executable.
	blastdbcmd_path : str
		Path to the blastdbcmd executable.

	Returns
	-------
	None
	"""
    
	input_paths = ff.read_lines(input_list)
	prefixes = get_file_prefixes(input_paths)
	# Create directory to store dummy data
	dummy_dir = os.path.join(output_directory, ct.DUMMY_DIR)
	ff.create_directory(dummy_dir)
	# Create dummy FASTA records
	dummy_seqids = [f'{i}-protein1' for i in (prefixes)]
	dummy_records = [ct.FASTA_RECORD_TEMPLATE.format(*[i, ct.DUMMY_PROT]) for i in dummy_seqids]
	dummy_fasta = os.path.join(dummy_dir, ct.DUMMY_FASTA)
	ff.write_lines(dummy_records, dummy_fasta)
	# Create BLAST db
	dummy_blastdb = os.path.join(dummy_dir, ct.DUMMY_BLASTDB)
	bf.make_blast_db(makeblastdb_path, dummy_fasta, dummy_blastdb, 'prot')
	# Get sequence from BLAST db
	dummy_blastdbcmd_fasta = os.path.join(dummy_dir, ct.DUMMY_BLASTDBCMD_FASTA)
	blastdbcmd_std = bf.run_blastdbcmd(blastdbcmd_path, dummy_blastdb, dummy_blastdbcmd_fasta)
	# Check if the sequence IDs in the BLAST db match the expected IDs
	records = SeqIO.parse(dummy_blastdbcmd_fasta, 'fasta')
	recids = [record.id for record in records]
	modified_seqids = set(dummy_seqids) - set(recids)
	# Delete dummy data
	shutil.rmtree(dummy_dir)
	# Exit if BLAST db includes sequence IDs that do not match the expected IDs
	if len(modified_seqids) > 0:
		pdb_prefix = [prefixes[p.split('-protein')[0]][0] for p in modified_seqids]
		sys.exit(ct.INPUTS_PDB_PREFIX.format('\n'.join(pdb_prefix)))


def adapt_loci(input_fastas: str, output_directory: str, cpu: int, bsr: float, translation_table: int) -> None:
    """

    Adapts an external schema for usage by calling chewBBACA PrepExternalSchema. 
    Removes invalid alleles and selects representative alleles to include in 
    the "short" directory.

    Parameters
    ----------
    input_fastas: str
        Path to the folder with the the fasta files. 
    output_directory: str
        Path to the directory where the final schema will be stored.
    cpu : int
        Number of CPU cores that will be used to run the process.
    blast_score_ratio : float
        The BLAST Score Ratio value that will be used to evaluate
        allele similarity and select representative alleles.
    translation_table : int
        Genetic code used to translate alleles.

    Returns
    -------
    None
        The function writes the output files to the specified directory.
    """

    pf.print_message("")
    pf.print_message("Starting External Schema Prep from chewBBACA...", "info")

    output_path = os.path.join(output_directory, "adapted_schema")

    cmd = [
        "chewBBACA.py",
        "PrepExternalSchema",
        "-g", input_fastas,
        "-o", output_path,
        "--cpu", str(cpu),
        "--bsr", str(bsr),
        "--t", str(translation_table)]

    # Run the PrepExternalSchema from chewie
    process = subprocess.Popen(
        cmd,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
        bufsize=1)

    # Create the output from chewBBACA to this log file
    for line in process.stdout:
        pf.print_message(line.strip(), "info")
    process.stdout.close()
    exit_code = process.wait()

    if exit_code == 0:
        pf.print_message("")
        pf.print_message("Schema creation completed", "info")
    else:
        pf.print_message("")
        pf.print_message(f"Schema creation failed with exit code {exit_code}", "error")

    # Validate schema fasta names
    pf.print_message('PREFIXES TRIAL START', 'info')
    gene_list_paths: str = os.path.join(output_directory, 'gene_list_paths.txt')
    with open(gene_list_paths, 'w') as gene_list:
        dir_list = os.listdir(input_fastas)
        for file in dir_list:
            if file.endswith(".fasta"):
                path = os.path.abspath(file)
                gene_list.write(f'{path}\n')
                
    makeblastdb_path = ff.join_paths(output_directory, [ct.MAKEBLASTDB_ALIAS])
    blastdbcmd_path = ff.join_paths(output_directory, [ct.BLASTDBCMD_ALIAS])
    pdb_prefixes = check_prefix_pdb(gene_list_paths, output_directory, makeblastdb_path, blastdbcmd_path)
    
    pf.print_message('PREFIXES TRIAL OVER', 'info')

