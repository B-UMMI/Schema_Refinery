#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Purpose
-------

Accepts a matrix with results from the AlleleCall process of
chewBBACA and determines the pairwise allelic differences to
create a distance matrix. It also determines the number of
shared loci to create a matrix with those values. The 'INF-'
prefix is removed and ASM, ALM, NIPH, NIPHEM, PLOT3, PLOT5,
LNF and LOTSC classifications are substituted by '0' before
determining pairwise distances.
"""


import os
import csv
import time
import shutil
import pickle
import random
import argparse
import traceback
import numpy as np
from multiprocessing import Pool

import mask_matrix as mm
import get_varSize_deep as gs


def pickle_dumper(content, output_file):
    """ Use the Pickle module to serialize an object.

        Parameters
        ----------
        content : type
            Variable that refers to the object that will
            be serialized and written to the output file.
        output_file : str
            Path to the output file.
    """

    with open(output_file, 'wb') as po:
        pickle.dump(content, po)


def pickle_loader(input_file):
    """ Use the Pickle module to de-serialize an object.

        Parameters
        ----------
        input_file : str
            Path to file with byte stream to be de-serialized.

        Returns
        -------
        content : type
            Variable that refers to the de-serialized
            object.
    """

    with open(input_file, 'rb') as pi:
        content = pickle.load(pi)

    return content


def progress_bar(process, total, tickval=5, ticknum=20, completed=False):
    """ Creates and prints progress bar to stdout.

        Parameters
        ----------
        process : multiprocessing.pool.MapResult
            Multiprocessing object.
        total : int
            Total number of inputs that have to be processed.
        tickval : int
            Progress completion percentage value for each
            tick.
        ticknum : int
            Total number of ticks in progress bar.
        completed : bool
            Boolean indicating if process has completed.

        Returns
        -------
        completed : bool
            Boolean indicating if process has completed.
    """

    # check if process has finished
    if (process.ready()):
        # print full progress bar and satisfy stopping condition
        progress_bar = '[{0}] 100%'.format('='*ticknum)
        completed = True

    # check how many inputs have been processed
    remaining = process._number_left
    if remaining == total:
        # print empty progress bar
        progress_bar = '[{0}] 0%'.format(' '*ticknum)
    else:
        # print progress bar, incremented by 5%
        progress = int(100-(remaining/total)*100)
        progress_tick = progress//tickval
        progress_bar = '[{0}{1}] {2}%'.format('='*progress_tick,
                                              ' '*(ticknum-progress_tick),
                                              progress)

    print('\r', progress_bar, end='')
    time.sleep(0.5)

    return completed


def function_helper(input_args):
    """ Runs function by passing set of provided inputs and
        captures exceptions raised during function execution.

        Parameters
        ----------
        input_args : list
            List with function inputs and function object to call
            in the last index.

        Returns
        -------
        results : list
            List with the results returned by the function.
            If an exception is raised it returns a list with
            the name of the function and the exception traceback.
    """

    try:
        results = input_args[-1](*input_args[0:-1])
    except Exception as e:
        func_name = (input_args[-1]).__name__
        traceback_lines = traceback.format_exception(etype=type(e), value=e,
                                                     tb=e.__traceback__)
        traceback_text = ''.join(traceback_lines)
        print('Error on {0}:\n{1}\n'.format(func_name, traceback_text))
        results = [func_name, traceback_text]

    return results


def map_async_parallelizer(inputs, function, cpu, callback='extend',
                           chunksize=1, show_progress=False):
    """ Parallelizes function calls by creating several processes
        and distributing inputs.

        Parameters
        ----------
        inputs : list
            List with inputs to process.
        function
            Function to be parallelized.
        cpu : int
            Number of processes to create (based on the
            number of cores).
        callback : str
            Results can be appended, 'append', to the
            list that stores results or the list of results
            can be extended, 'extend'.
        chunksize : int
            Size of input chunks that will be passed to
            each process. The function will create groups
            of inputs with this number of elements.
        show_progress: bool
            True to show a progress bar with the percentage
            of inputs that have been processed, False
            otherwise.

        Returns
        -------
        results : list
            List with the results returned for each function
            call.
    """

    results = []
    pool = Pool(cpu)
    if callback == 'extend':
        rawr = pool.map_async(function, inputs,
                              callback=results.extend, chunksize=chunksize)
    elif callback == 'append':
        rawr = pool.map_async(function, inputs,
                              callback=results.append, chunksize=chunksize)

    if show_progress is True:
        completed = False
        while completed is False:
            completed = progress_bar(rawr, len(inputs))

    rawr.wait()

    return results


def determine_distances_fast_parallel(input_table, cpu_cores, tmp_directory, genome_ids):
    """
    """

    # import matrix without column and row identifiers
    with open(input_table, 'r') as infile:
        lines = ('\t'.join(line.split('\t')[1:]) for line in infile if not line.startswith('FILE'))
        # dtype=float32 should be faster than integer dtypes
        # but runs faster with dtype=int32 in test setup
        # dtype=int32 supports max integer of 2147483647
        # should be safe even when arrays are multiplied
        np_matrix = np.loadtxt(fname=lines, delimiter='\t', dtype='int32')

    rows_indexes = [i for i in range(len(np_matrix))]
    random.shuffle(rows_indexes)
    # divide inputs into 20 lists for 5% progress resolution
    parallel_inputs = divide_list_into_n_chunks(rows_indexes, 20)

    common_args = [[l, np_matrix, genome_ids, tmp_directory, separate_dists] for l in parallel_inputs]

    # increasing cpu cores can greatly increase memory usage
    results = map_async_parallelizer(common_args,
                                     function_helper,
                                     cpu_cores,
                                     show_progress=True)

    return results


def separate_dists(indexes, np_matrix, genome_ids, tmp_directory):

    # multiply one row per cycle to avoid memory overflow
    output_files = {}
    for i in indexes:
        current_genome = genome_ids[i]
        # get one row to perform pairwise comparisons against whole matrix
        current_row = np_matrix[i:i+1, :]
        # do not multiply against rows that were multiplied against
        # matrix's rows in previous iterations
        # combinations instead of permutations
        permutation_rows = np_matrix[i:, :]

        # multiply 1D-array per whole matrix
        # all non-shared loci will be converted to 0
        # values different than 0 correspond to shared loci
        multiplied = current_row * permutation_rows
        # count number of shared loci, non-zero values
        pairwise_shared_loci = np.count_nonzero(multiplied, axis=-1)

        # subtraction will lead to values different than 0 for loci that have different alleles
        # multiplying ensures that we only keep results for shared loci and not for
        # loci that are not shared and that had value different than 0 from subtraction
        pairwise_allelic_differences = np.count_nonzero(multiplied * (current_row - permutation_rows), axis=-1)

        output_file = os.path.join(tmp_directory, current_genome)
        pickle_dumper([pairwise_shared_loci, pairwise_allelic_differences], output_file)
        output_files[current_genome] = output_file

    return output_files


def divide_list_into_n_chunks(list_to_divide, n):
    """ Divides a list into a defined number of sublists.

        Parameters
        ----------
        list_to_divide : list
            List to divide into sublists.
        n : int
            Number of sublists to create.

        Returns
        -------
        sublists : list
            List with the sublists created by dividing
            the input list.
    """

    sublists = []
    d, r = divmod(len(list_to_divide), n)
    for i in range(n):
        si = (d+1)*(i if i < r else r) + d*(0 if i < r else i - r)
        sublists.append(list_to_divide[si:si+(d+1 if i < r else d)])

    # exclude lists that are empty due to small number of elements
    sublists = [i for i in sublists if len(i) > 0]

    return sublists
    

def join_iterable(iterable, delimiter='\t'):
    """
    """

    joined = delimiter.join(iterable)

    return joined


def write_text(text, output_file, mode='w'):
    """
    """

    # write matrix to output file
    with open(output_file, mode) as outfile:
        outfile.write(text+'\n')


def write_lines(lines, output_file, mode='w'):
    """ Writes a matrix to a file.

    Parameters
    ----------
    matrix_rows : list
        List of sublists where each sublist corresponds
        to one row of a AlleleCall matrix.
    output_matrix : str
        Path to the file that should be created to store
        the output matrix.

    Returns
    -------
    Writes matrix rows (each sublist of the input list is
    a row) to the output file.
    """

    # join matrix lines into chunk of text
    concat_lines = [join_iterable(line, '\t')
                    for line in lines]
    lines_text = join_iterable(concat_lines, '\n')

    write_text(lines_text, output_file, mode)


def read_tabular(input_file, delimiter='\t'):
    """
    Read a tabular file.

    Parameters
    ----------
    input_file : str
        Path to a tabular file.
    delimiter : str
        Delimiter used to separate file fields.

    Returns
    -------
    lines : list
        A list with a sublist per line in the input file.
        Each sublist has the fields that were separated by
        the defined delimiter.
    """

    with open(input_file, 'r') as infile:
        reader = csv.reader(infile, delimiter=delimiter)
        lines = [line for line in reader]

    return lines


def get_input_ids(input_file, delimiter='\t'):

    with open(input_file, 'r') as infile:
        reader = csv.reader(infile, delimiter=delimiter)
        input_ids = [line[0] for line in reader][1:]

    return input_ids


def merge_dictionaries(dictionaries):
    """
    """

    merged = {}
    for d in dictionaries:
        merged = {**merged, **d}

    return merged


def main(input_matrix, output_directory, cpu_cores):

    # create output directory if it does not exist
    if os.path.isdir(output_directory) is False:
        os.mkdir(output_directory)

    # determine input basename
    input_basename = os.path.basename(input_matrix)
    # remove extension that is after last '.'
    input_basename = '.'.join(input_basename.split('.')[0:-1])

    # define '0' as masking characters for all non-numeric
    # classifications
    classes = ['ALM', 'ASM', 'LNF', 'NIPH',
               'NIPHEM', 'PLOT3', 'PLOT5', 'LOTSC']
    masking_dict = {c: '0' for c in classes}

    print('Masking matrix before determining pairwise distances...', end='')
    output_masked = os.path.join(output_directory,
                                 '{0}_masked.tsv'.format(input_basename))
    total_masked = mm.mask_matrix(input_matrix, masking_dict, output_masked)
    # mask matrix
    print('masked matrix available at {0}'.format(output_masked))

    # create temp directory to store pairwise distances per genome
    tmp_directory = os.path.join(output_directory, 'tmp')
    if os.path.isdir(tmp_directory) is False:
        os.mkdir(tmp_directory)

    # get sample identifiers
    genome_ids = get_input_ids(input_matrix, delimiter='\t')
    total_genomes = len(genome_ids)

    # this uses a lot of memory for really large matrices
    print('Determining pairwise distances...')
    start2 = time.time()
    parallel_pairwise = determine_distances_fast_parallel(output_masked, cpu_cores, tmp_directory, genome_ids)
    end2 = time.time()
    delta2 = end2 - start2
    print('\n', delta2)

    merged = merge_dictionaries(parallel_pairwise)

    print('Creating distance matrix...', end='')
    # create files with headers
    headers = [''] + genome_ids
    headers = join_iterable(headers, '\t')
    output_pairwise = os.path.join(output_directory,
                                   '{0}_allelic_differences.tsv'.format(input_basename))
    write_text(headers, output_pairwise)

    output_p = os.path.join(output_directory,
                            '{0}_shared_loci.tsv'.format(input_basename))
    write_text(headers, output_p)

    # import arrays per genome and save to matrix file
    sl_lines = []
    ad_lines = []
    limit = 500
    for g in genome_ids:
        current_file = merged[g]
        # load data
        data = pickle_loader(current_file)

        shared_loci = list(data[0])
        shared_loci = list(map(str, shared_loci))
        allele_diffs = list(data[1])
        allele_diffs = list(map(str, allele_diffs))

        padding = [''] * (len(genome_ids)-len(allele_diffs))

        sl_line = [g] + padding + shared_loci
        sl_lines.append(sl_line)
        ad_line = [g] + padding + allele_diffs
        ad_lines.append(ad_line)

        if len(sl_lines) >= limit or g == genome_ids[-1]:
            write_lines(ad_lines, output_pairwise, mode='a')
            ad_lines = []
            write_lines(sl_lines, output_p, mode='a')
            sl_lines = []

    print('done.')
    print('Distance matrix with allelic differences available at {0}'.format(output_pairwise))
    print('Matrix with shared loci available at {0}'.format(output_p))

    # delete folder with intermediate pickles
    shutil.rmtree(tmp_directory)


def parse_arguments():

    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument('-i', '--input-matrix', type=str,
                        required=True, dest='input_matrix',
                        help='Path to the file with the AlleleCall '
                             'matrix (default name given by chewBBACA'
                             ' is results_alleles.tsv).')

    parser.add_argument('-o', '--output-directory', type=str,
                        required=True, dest='output_directory',
                        help='Path to the output directory that will be '
                             'created to store the results.')

    parser.add_argument('-c', '--cpu-cores', type=int,
                        required=False, default=1,
                        dest='cpu_cores',
                        help='Number of CPU cores used to perform '
                             'pairwise comparisons.')

    args = parser.parse_args()

    return args


if __name__ == "__main__":

    args = parse_arguments()
    main(**vars(args))
