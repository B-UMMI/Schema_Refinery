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


try:
    from utils import (print_functions as pf)
except ModuleNotFoundError:
    from SchemaRefinery.utils import (print_functions as pf)

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

    cmd = [
        "chewBBACA.py",
        "PrepExternalSchema",
        "-g", input_fastas,
        "-o", output_directory,
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

