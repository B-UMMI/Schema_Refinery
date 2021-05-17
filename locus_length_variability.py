#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Purpose
-------


Code documentation
------------------
"""


import os
import argparse
import statistics
from collections import Counter

from Bio import SeqIO
import plotly.graph_objs as go
from plotly.offline import plot
from plotly.subplots import make_subplots


def import_sequences(fasta_path):
    """ Imports sequences from a FASTA file.

        Parameters
        ----------
        fasta_path : str
            Path to a FASTA file.

        Returns
        -------
        seqs_dict : dict
            Dictionary that has sequences ids as keys and
            sequences as values.
    """

    records = SeqIO.parse(fasta_path, 'fasta')
    seqs_dict = {rec.id: str(rec.seq.upper()) for rec in records}

    return seqs_dict


def histogram_tracer(xdata, ydata, marker_color):
    """
    """

    trace = go.Bar(x=xdata,
                   y=ydata,
                   marker_color=marker_color)

    return trace


def create_subplots(nplots, ncols, titles):
    """
    """

    cols = ncols
    rows = statistics.math.ceil((nplots / ncols))
    subplots_fig = make_subplots(rows=rows, cols=cols,
                                 subplot_titles=titles)

    return subplots_fig


def add_multiple_traces(traces, fig, srow, scol, max_cols):
    """
    """

    row_start = srow
    col_start = scol
    for k, v in traces.items():
        trace = v[1]
        fig.append_trace(trace, row_start, col_start)
        if col_start < max_cols:
            col_start += 1
        else:
            col_start = 1
            row_start += 1

    return fig


def write_text(text, output_file, end='\n'):
    """
    """

    with open(output_file, 'w') as outfile:
        outfile.write(text+end)


def write_stats(stats, header, output_file):
    """
    """

    lines = [header]
    for locus, values in stats.items():
        locus_id = (os.path.basename(locus)).split('.fasta')[0]
        most_frequent_str = ';'.join([str(v) for v in values[-1]])
        str_values = [str(v) for v in values[0:-1]] + [most_frequent_str]
        current = [locus_id] + str_values
        lines.append('\t'.join(current))

    text = '\n'.join(lines)
    write_text(text, output_file)


def main(schema_directory, output_directory):

    if os.path.isdir(output_directory) is not True:
        os.mkdir(output_directory)

    # list genes in schema
    schema_loci = [file
                   for file in os.listdir(schema_directory)
                   if '.fasta' in file]
    schema_loci = [os.path.join(schema_directory, file)
                   for file in schema_loci]

    full_stats = {}
    single_mode_ids = []
    single_mode_stats = {}
    single_mode_tracer = {}
    multi_mode_ids = []
    multi_mode_stats = {}
    multi_mode_tracer = {}
    for locus in schema_loci:
        records = import_sequences(locus)

        # get total number of records
        total_alleles = len(records)

        # get length of each allele
        allele_lengths = {recid: len(seq) for recid, seq in records.items()}
        lengths_list = list(allele_lengths.values())

        # determine mean
        locus_mean = round(statistics.mean(lengths_list), 1)

        # determine median
        locus_median = round(statistics.median(lengths_list), 1)

        # determine mode and create histograms for loci with single mode
        try:
            locus_mode = statistics.mode(lengths_list)
        except Exception as e:
            print(locus, e)
            locus_mode = 'None'

        # determine if there are other values that are very frequent
        length_counts = Counter(lengths_list)
        # get values ordered starting with most frequent
        freq_ordered = length_counts.most_common()
        most_frequent = freq_ordered[0]
        # select length values that are also very frequent
        high_freq = [v
                     for v in freq_ordered[1:]
                     if v[1] >= (0.3*most_frequent[1])]
        frequent_values = [most_frequent] + high_freq

        # select length values that are also very frequent
        asm = most_frequent[0]-(most_frequent[0]*0.2)
        alm = most_frequent[0]+(most_frequent[0]*0.2)
        outlier = [v
                   for v in freq_ordered[1:]
                   if (v[0] >= alm or v[0] <= asm)]

        xdata = [v[0] for v in freq_ordered]
        ydata = [v[1] for v in freq_ordered]

        if locus_mode != 'None' and len(outlier) > 0:
            # create histogram object
            trace = histogram_tracer(xdata, ydata, '#0868ac')
            single_mode_tracer[locus] = [frequent_values, trace]
            single_mode_ids.append((os.path.basename(locus)).split('.fasta')[0])
            # add data to dict
            single_mode_stats[locus] = [total_alleles, locus_mean,
                                        locus_median, locus_mode,
                                        frequent_values]
        elif locus_mode == 'None':
            trace = histogram_tracer(xdata, ydata, '#0868ac')
            multi_mode_tracer[locus] = [frequent_values, trace]
            multi_mode_ids.append((os.path.basename(locus)).split('.fasta')[0])
            # add data to dict
            multi_mode_stats[locus] = [total_alleles, locus_mean, locus_median,
                                       locus_mode, frequent_values]

        full_stats[locus] = [total_alleles, locus_mean,
                             locus_median, locus_mode,
                             frequent_values]

    # header for file with statistics
    header = 'locus\tnum_alleles\tmean\tmedian\tmode\tmost_frequent'
    # template for file path
    file_template = '{0}/{1}'

    # save full stats
    full_stats_outfile = file_template.format(output_directory,
                                              'full_stats.tsv')
    write_stats(full_stats, header, full_stats_outfile)

    # create TSV lines with single mode loci data
    single_mode_outfile = file_template.format(output_directory,
                                               'single_mode_stats.tsv')

    write_stats(single_mode_stats, header, single_mode_outfile)

    # write single mode high length variability identifiers to file
    single_mode_ids_outfile = file_template.format(output_directory,
                                                   'single_mode_ids.txt')
    single_mode_ids_text = '\n'.join(single_mode_ids)
    write_text(single_mode_ids_text, single_mode_ids_outfile)

    # create HTML with subplots for loci with single mode
    titles = [(os.path.basename(k)).split('.fasta')[0] for k in single_mode_tracer]
    single_mode_fig = create_subplots(len(single_mode_tracer), 2, titles)
    single_mode_fig = add_multiple_traces(single_mode_tracer, single_mode_fig,
                                          1, 1, 2)

    fig_height = (statistics.math.ceil(len(single_mode_fig.data)/2))*300
    single_mode_fig.update_layout(title='Loci with alleles 20% smaller or larger than the allele length mode',
                                  height=fig_height,
                                  showlegend=False,
                                  bargap=0.1,
                                  template='ggplot2')

    single_mode_fig.update_xaxes(type='log', title_text='Allele length (log)')
    single_mode_fig.update_yaxes(type='log', title_text='Number of alleles (log)')

    single_mode_plotfile = file_template.format(output_directory,
                                                'single_mode.html')
    plot(single_mode_fig, filename=single_mode_plotfile, auto_open=True)

    # create TSV lines with multi mode loci data
    multi_mode_outfile = file_template.format(output_directory,
                                              'multi_mode_stats.tsv')
    write_stats(multi_mode_stats, header, multi_mode_outfile)

    # write multi mode identifiers to file
    multi_mode_ids_outfile = file_template.format(output_directory,
                                                  'multi_mode_ids.txt')
    multi_mode_ids_text = '\n'.join(multi_mode_ids)
    write_text(multi_mode_ids_text, multi_mode_ids_outfile)

    # create HTML with subplots for loci with multiple modes
    titles = [(os.path.basename(k)).split('.fasta')[0] for k in multi_mode_tracer]
    multimodal_fig = create_subplots(len(multi_mode_tracer), 2, titles)
    multimodal_fig = add_multiple_traces(multi_mode_tracer, multimodal_fig,
                                         1, 1, 2)

    fig_height = (statistics.math.ceil(len(multimodal_fig.data)/2))*300
    multimodal_fig.update_layout(title='Loci with multimodal allele length distributions',
                                 height=fig_height,
                                 showlegend=False,
                                 bargap=0.1,
                                 template='ggplot2')

    multimodal_fig.update_xaxes(type='log', title_text='Allele length (log)')
    multimodal_fig.update_yaxes(type='log', title_text='Number of alleles (log)')

    plot_file = file_template.format(output_directory,
                                     'multiple_mode.html')
    plot(multimodal_fig, filename=plot_file, auto_open=True)


def parse_arguments():

    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument('-s', type=str, required=True,
                        dest='schema_directory',
                        help='')

    parser.add_argument('-o', type=str, required=True,
                        dest='output_directory',
                        help='')

    args = parser.parse_args()

    return args


if __name__ == '__main__':

    args = parse_arguments()
    main(**vars(args))
