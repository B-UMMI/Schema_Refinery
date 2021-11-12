#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
AUTHOR

    Pedro Cerqueira
    github: @pedrorvc

DESCRIPTION

    This script serves to analyse assemblies generated from a pipeline like Innuca 
    and create a report based on that analysis.


    USAGE:

    # Report generation
    post_innuca_report.py -i path/to/assemblies/dir
                            -o path/to/output/dir --cpu 6 --total_bps 2300000
                            --nr_contigs 350 --gc_content 0.45
"""

import os
import bz2
import gzip
import time
import zipfile
import argparse
import itertools
import statistics as stats
from multiprocessing import Pool, cpu_count

import pandas as pd

from Bio import SeqIO
from Bio.SeqUtils import GC


COPEN = {
    "gz": gzip.open,
    "bz2": bz2.open,
    "zip": zipfile.ZipFile
}

MAGIC_DICT = {
    b"\x1f\x8b\x08": "gz",
    b"\x42\x5a\x68": "bz2",
    b"\x50\x4b\x03\x04": "zip"
}

"""
dict: Dictionary containing the binary signatures for three compression formats
(gzip, bzip2 and zip).
"""


def guess_file_compression(file_path, magic_dict=None):
    """Guesses the compression of an input file.
    This function guesses the compression of a given file by checking for
    a binary signature at the beginning of the file. These signatures are
    stored in the :py:data:`MAGIC_DICT` dictionary. The supported compression
    formats are gzip, bzip2 and zip. If none of the signatures in this
    dictionary are found at the beginning of the file, it returns ``None``.
    
    Args:
        file_path : str
            Path to input file.
        magic_dict : dict, optional
            Dictionary containing the signatures of the compression types. The
            key should be the binary signature and the value should be the
            compression format. If left ``None``, it falls back to
            :py:data:`MAGIC_DICT`.
    Returns:
        file_type : str or None
            If a compression type is detected, returns a string with the format.
            If not, returns ``None``.
    """

    if not magic_dict:
        magic_dict = MAGIC_DICT

    max_len = max(len(x) for x in magic_dict)

    with open(file_path, "rb") as f:
        file_start = f.read(max_len)

    for magic, file_type in magic_dict.items():
        if file_start.startswith(magic):
            return file_type

    return None


def is_fasta(filename):
    """ Checks if a file is a FASTA file.

        Args:
            filename (str): the full path to the FASTA file

        Returns:
            True if FASTA file,
            False otherwise
    """
    
    
    try:
        with open(filename, "r") as handle:
            fasta = SeqIO.parse(handle, "fasta")
    
            # returns True if FASTA file, False otherwise
            return any(fasta)
    
    except UnicodeDecodeError:
        return False


def is_fasta_gz(filename):
    """ Checks if a file is a FASTA GZ file.
    
        Args:
            filename (str): the full path to the FASTA file

        Returns:
            True if FASTA file,
            False otherwise
    """
    
    with gzip.open(filename, "rt") as handle:
        fasta_gz = SeqIO.parse(handle, "fasta")
    
        # returns True if FASTA file, False otherwise
        return any(fasta_gz)


def check_if_list_or_folder(folder_or_list):
    """ Checks if the input is a file or a directory.

        Args:
            folder_or_list (str): the full path to the file or directory

        Returns:
            list_files (str) if folder_or_list is a path to a file,
            list_files (list) if folder_or_list is a path to a directory,
            Raises Exception otherwise
    """

    # check if input argument is a file or a directory
    if os.path.isfile(folder_or_list):
        list_files = folder_or_list

    elif os.path.isdir(folder_or_list):

        fasta_files = []

        for genome in os.listdir(folder_or_list):

            genepath = os.path.join(folder_or_list, genome)

            # do not continue if genepath is a dir
            if os.path.isdir(genepath):
                continue

            # check if file is a FASTA file
            if is_fasta(genepath):
                fasta_files.append(os.path.abspath(genepath))
            
            # check if file is a FASTA file
            elif is_fasta_gz(genepath):
                fasta_files.append(os.path.abspath(genepath))

        # if there are FASTA files
        if fasta_files:
            # store full paths to FASTA files
            with open("listGenes.txt", "w") as f:
                for genome in fasta_files:
                    f.write(genome + "\n")
        else:
            raise Exception("There were no FASTA files in the given directory. Please provide a directory \
                            with FASTA files or a file with the list of full paths to the FASTA files.")

        list_files = "listGenes.txt"

    else:
        raise Exception("Input argument is not a valid directory or file with a list of paths. \
                        Please provide a valid input, either a folder with FASTA files or a file with \
                        the list of full paths to FASTA files (one per line).")

    return list_files


def verify_cpu_usage(cpu_to_use):
    """ Verify the cpu usage for the script.

        Args:
            cpu_to_use (int): the number of cpu provided to post_innuca_report

        Returns:
            cpu_to_use (int): the number of cpu to use after verification

        Example:

            >>> verify_cpu_usage(6)
            6
    """
    total_cpu = cpu_count()

    # do not allow a value of cpuToUse greater than the number of cores/threads
    if cpu_to_use > total_cpu:
        print("Warning! You have provided a CPU core count value that exceeds the number of cores in your machine!")
        print("Setting a different value for the CPU core count...")
        # define a value that is safe according to the number of available cores/threads
        if total_cpu > 2:
            cpu_to_use = total_cpu - 2
        elif total_cpu == 2:
            cpu_to_use = 1
        print("CPU core count value set to: ", cpu_to_use)

    elif cpu_to_use < total_cpu and cpu_to_use > total_cpu - 2:
        print("Warning! You have provided a CPU core count value that is close to the maximum core count of your machine ("
              + str(cpu_to_use) + '/' + str(total_cpu) + "). This may affect your system responsiveness.")

    return cpu_to_use


def track_job(job, update_interval=3):
    """ Tracks multiprocessing jobs
    """
    
    while job._number_left > 0:
        print("Tasks remaining = {0}".format(
        job._number_left * job._chunksize))
        time.sleep(update_interval)


def flatten_list(list_to_flatten):
    """Flattens one level of a nested list

        Args:
            list_to_flatten (list)

        Returns:
            flattened list

        Example:

            >>> flatten_list([[[1,2],[3,4]]])
            [[1, 2], [3, 4]]

    """

    return list(itertools.chain(*list_to_flatten))


def analyse_report(report, nr_contigs, min_bps, max_bps, min_gc, max_gc):
    """ Classify samples using report results

        Args:
            report (dict): Contains the assembly and annotation reports
            total_bps (int): Total assembled base pairs
            nr_contigs (int): number of generated contigs
            gc_content (float): Percentage of GC content allowed

        Returns:

            result (dict) : Contains the results of the analysis
    """

    result = {}

    for record in report:

        result[record["Sample"]] = []

        if record["Total assembly length"] < min_bps:
            result[record["Sample"]].append(["FAIL", "Low_BPs"])
        elif record["Total assembly length"] > max_bps:
            result[record["Sample"]].append(["FAIL", "High_BPs"])

        if record["GC content"] < min_gc:
            result[record["Sample"]].append(["FAIL", "Low_GC"])
        elif record["GC content"] > max_gc:
            result[record["Sample"]].append(["FAIL", "High_GC"])

        if record["Number of contigs"] > nr_contigs:
            result[record["Sample"]].append(["FAIL", "Nr_contigs"])

        if result[record["Sample"]] == []:
            result[record["Sample"]].append("PASS")

    return result


def write_final_report(report, analysed_report, output):
    """ Writes a tab-delimited file of analysis results.

        Args:
            report (pandas.DataFrame): DataFrame containing the assembly and annotation results
            analysed_report (dict): Contains the results of the report analysis
            output (str): Path to the output directory

        Returns:
            Tab-delimited file in the specified output director
    """

    # Build the final report (DataFrame)
    final_report = report.copy()

    # Convert analysed report (dict) to pandas.Series
    s = pd.Series(analysed_report, name='Result')

    # Change the name of the index
    s.index.name = 'Sample'

    # Reset index
    s = s.reset_index()

    # Add the result Series to the final_report DataFrame
    final_report["Result"] = s["Result"]

    # Write file
    final_report.to_csv(os.path.join(output, "final_report.tsv"),
                        sep='\t',
                        encoding='utf-8',
                        index=False)

    print("File written to {}".format(output))


def calc_n50(contig_sizes):
    """ Calculates the N50 of an assembly.

        Args:
            contig_sizes (list): Contains the sizes of the contigs from an assembly file

        Returns:
            l (int): Calculated n50
    """
    
    # Sort the contig sizes in descending order
    contig_sizes.sort(reverse=True)
    
    # Calculate n50
    s = sum(contig_sizes)
    limit = s * 0.5
    for l in contig_sizes:
        s -= l
        if s <= limit:
            return l

        
def analyse_assembly(assembly):
    """ Analyses an assembly file

        Args:
            assemblies (str): Path to the assembly file

        Returns:
            results (dict): Contains the results of the analysis
    """
    
    
    assembly_file = assembly
  
    # Get the sample name from the file
    sample = os.path.basename(assembly_file).split(".")[0]
    
    # Guess the file compression
    ftype = guess_file_compression(assembly_file)
    
    # File is not compressed
    if ftype is None:    
        # Get the records of the assembly file
        records = list(SeqIO.parse(assembly_file, "fasta"))
    
    # This can guess the compression of gz, bz2 and zip.
    else:      
        records = []
        
        with COPEN[ftype](assembly_file, "rt") as af:
            for record in SeqIO.parse(af, "fasta"):
                records.append(record)
    
    # Calculate the GC content
    all_gc_content = [GC(seq.seq) for seq in records]
    gc_content = stats.mean(all_gc_content) / 100

    # Get the total number of contigs in the assembly file
    nr_contigs = len(records)

    # Get the contig sizes
    sizes = [len(seq) for seq in records]

    # Calculate the average contig size
    avg_size = stats.mean(sizes)

    # Calculate the total assembly length
    total_length = sum(sizes)

    # Calculate the N50
    n50 = calc_n50(sizes)

    # Determine missing data
    proportion_missing_data = (sum([rec.seq.count("N") for rec in records]) / 
                               total_length)

    # Save the results in a dictionary
    results = {"Sample": sample, "Number of contigs": nr_contigs,
               "Average contig size": avg_size, "N50": n50,
               "Total assembly length": total_length, "GC content": gc_content,
               "Missing Data": round(proportion_missing_data, ndigits=4)}

    return results


def main(output, assemblies, cpu, nr_contigs, min_bps, max_bps, min_gc, max_gc):

    # avoid user to run the script with all cores available,
    # could impossibilitate any usage when running on a laptop
    cpu_to_apply = verify_cpu_usage(cpu)

    # check if output directory exists
    if not os.path.exists(output):
        os.mkdir(output)

    if assemblies:

        assemblies_file = check_if_list_or_folder(assemblies)


        listGenes = []
        with open(assemblies_file, "r") as gf:
            for gene in gf:
                gene = gene.rstrip('\n')
                listGenes.append(gene)
        listGenes.sort()

        try:
            os.remove("listGenes.txt")
        except Exception:
            pass
        
        # List to save the results of the multiprocessing
        assembly_analysis_results = []

        print("Calculating assembly statistics...\n")
        p = Pool(processes=cpu_to_apply)
        r = p.map_async(analyse_assembly, listGenes, callback=assembly_analysis_results.append)
        
        track_job(r, 3)

        r.wait()

        print("\nAnalysing results...\n")

        # Flatten result nested list by one level
        results = flatten_list(assembly_analysis_results)

        # Analyse results
        analysed_report = analyse_report(results, nr_contigs,
                                         min_bps, max_bps, min_gc, max_gc)

        # Print amount of Fails to console
        fails = sum([fail[0].count("FAIL") for fail in analysed_report.values()])

        print("The analysis detected {0} FAILS on a total of {1} assemblies. Check the report for more details.\n".format(fails, len(listGenes)))

        print("Writing report...\n")

        # Convert dictionary into pandas DataFrame
        report = pd.DataFrame(results)
        report = report[['Sample', 'Number of contigs', 'Average contig size', 'N50',
                         'Total assembly length', 'GC content', 'Missing Data']]

        # Write the final report
        write_final_report(report, analysed_report, output)

    print("Execution Finished")


def parse_arguments():

    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument('-o', '--output', type=str, required=True, dest='output_path',
                        default=False, help='Path to the output directory.')

    parser.add_argument('-i', '--assemblies', type=str, required=False, dest='assembly_path',
                        default=False, help='Path to the directory containing the assemblies.')


    parser.add_argument('--cpu', type=int, required=False, dest='cpu',
                        default=2, help='Number of CPU cores to use. If the provided value exceeds \
                                         the maximum number of available cores uses maximum -2.')

    parser.add_argument('--nr_contigs', type=int, required=False, dest='nr_contigs',
                        default=350, help='Number of contigs allowed for each assembly.')

    parser.add_argument('--min_bps', type=int, required=False, dest='minimum_number_of_bases',
                        default=1, help='Minimum number of total bases accepted for a genome/assembly.')

    parser.add_argument('--max_bps', type=int, required=False, dest='maximum_number_of_bases',
                        default=9999999999999999, help='Maximum number of total bases accepted for a genome/assembly.')

    parser.add_argument('--min_gc', type=float, required=False, dest='minimum_gc_content',
                        default=0.0, help='Minimum GC content value.')

    parser.add_argument('--max_gc', type=float, required=False, dest='maximum_gc_content',
                        default=1.0, help='Minimum GC content value.')


    args = parser.parse_args()

    return [args.output_path, args.assembly_path, args.cpu, args.nr_contigs,
            args.minimum_number_of_bases, args.maximum_number_of_bases,
            args.minimum_gc_content, args.maximum_gc_content]


if __name__ == '__main__':

    args = parse_arguments()
    main(args[0], args[1], args[2], args[3], args[4], args[5], args[6], args[7])
