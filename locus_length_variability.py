#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Purpose
-------
This scripts computes allele length statistics for all
loci in a schema, identifying loci with multimodal allele
length distributions and loci with alleles that deviate
from the mode value.

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
    """Import sequences from a FASTA file.

    Parameters
    ----------
    fasta_path : str
        Path to a FASTA file.

    Returns
    -------
    seqs_dict : dict
        Dictionary that has sequence identifiers as
        keys and sequences as values.
    """
    records = SeqIO.parse(fasta_path, 'fasta')
    seqs_dict = {rec.id: str(rec.seq.upper()) for rec in records}

    return seqs_dict


def histogram_trace(xdata, ydata, marker_colors):
    """Create a histogram (go.Bar) trace.

    Parameters
    ----------
    xdata : list
        Sets the x coordinates (sequence
        length values).
    ydata : list
        Sets the y coordinates (counts for
        each sequence length value).
    marker_colors : str
        Colors to fill the bars.

    Returns
    -------
    trace : plotly.graph_objects.Bar
        Trace with data to display histogram.
    """
    trace = go.Bar(x=xdata,
                   y=ydata,
                   marker_color=marker_colors)

    return trace


def create_subplots(nplots, ncols, titles):
    """Create a Figure objcet with predefined subplots.

    Parameters
    ----------
    nplots : int
        Number of plots that will be displayed.
    ncols : int
        Number of subplots per row.
    titles : list
        Subplot titles.

    Returns
    -------
    subplots_fig : plotly.graph_objs._figure.Figure
        Figure object with predefined subplots.
    """
    # determine number of rows
    nrows = statistics.math.ceil((nplots / ncols))
    subplots_fig = make_subplots(rows=nrows, cols=ncols,
                                 subplot_titles=titles)

    return subplots_fig


def add_multiple_traces(traces, fig, max_cols):
    """Add traces to Figure object based on the number of cols per row.

    Parameters
    ----------
    traces : list
        List with trace objects to add.
    fig : plotly.graph_objs._figure.Figure
        Figure object with predefined subplots.
    max_cols : int
        Maximum number of columns per row.

    Returns
    -------
    fig : plotly.graph_objs._figure.Figure
        Figure object with the data to display the
        traces.
    """
    row_start = 1
    col_start = 1
    for t in traces:
        fig.append_trace(t, row_start, col_start)
        if col_start < max_cols:
            col_start += 1
        else:
            col_start = 1
            row_start += 1

    return fig


def write_text(text, output_file, end='\n'):
    """Write text to a file.

    Parameters
    ----------
    text : str
        Text to write.
    output_file : str
        Path to output file.
    end : str
        Characters to add to end of file.
    """
    with open(output_file, 'w') as outfile:
        outfile.write(text+end)


def main(schema_directory, output_directory, size_threshold, frequency_threshold):

    if os.path.isdir(output_directory) is not True:
        os.mkdir(output_directory)

    # list genes in schema
    schema_loci = [file
                   for file in os.listdir(schema_directory)
                   if '.fasta' in file]
    schema_loci = [os.path.join(schema_directory, file)
                   for file in schema_loci]

    loci_stats = {}
    for locus in schema_loci:
        records = import_sequences(locus)

        # get total number of records
        total_alleles = len(records)

        # get length of each allele
        allele_lengths = {recid: len(seq) for recid, seq in records.items()}
        lengths_list = list(allele_lengths.values())

        locus_id = os.path.basename(locus).split('.fasta')[0]

        # determine sequence length mean
        locus_mean = round(statistics.mean(lengths_list), 1)

        # determine sequence length median
        locus_median = round(statistics.median(lengths_list), 1)

        # determine mode and create histograms for loci with single mode
        try:
            locus_mode = statistics.mode(lengths_list)
        # locus has at least two length values that are equally frequent
        except Exception as e:
            locus_mode = 'None'

        # determine if there are other values that are very frequent
        length_counts = Counter(lengths_list)
        # get values ordered starting with most frequent
        freq_ordered = length_counts.most_common()
        most_frequent = freq_ordered[0]
        # select length values that are also very frequent (>=0.3)
        high_freq = [v
                     for v in freq_ordered[1:]
                     if v[1] >= (frequency_threshold*most_frequent[1])]

        frequent_values = [most_frequent] + high_freq

        # determine if there are alleles outside (mode +/- (mode*size_threshold))
        asm = most_frequent[0]-(most_frequent[0]*size_threshold)
        alm = most_frequent[0]+(most_frequent[0]*size_threshold)
        outlier = [v
                   for v in freq_ordered[1:]
                   if (v[0] >= alm or v[0] <= asm)]

        # sort by length value
        sorted_freq = sorted(freq_ordered, key=lambda x: x[0])
        xdata = [v[0] for v in sorted_freq]
        ydata = [v[1] for v in sorted_freq]
        colors = []
        for v in xdata:
            if v <= asm:
                colors.append('#fec44f')
            elif v >= alm:
                colors.append('#bd0026')
            elif v == most_frequent[0]:
                colors.append('#238b45')
            elif len(high_freq) > 0 and v in [f[0] for f in high_freq]:
                colors.append('#0868ac')
            else:
                colors.append('#969696')

        loci_stats[locus] = [total_alleles, locus_mean,
                             locus_median, locus_mode,
                             frequent_values, locus_id]

        # create histogram for loci that have multimodal distribution or
        # loci that have alleles outside the acceptable length interval
        if locus_mode == 'None' or len(outlier) > 0:
            # create histogram object
            trace = histogram_trace(xdata, ydata, colors)
            if len(outlier) > 0 and locus_mode == 'None':
                loci_stats[locus].extend([trace, 'outlier+multimodal'])
            elif len(outlier) > 0:
                loci_stats[locus].extend([trace, 'outlier'])
            elif locus_mode == 'None':
                loci_stats[locus].extend([trace, 'multimodal'])
        else:
            loci_stats[locus].extend([None, 'single'])

    # header for file with statistics
    header = ('locus\tnum_alleles\tmean\tmedian\tmode\t'
              'most_frequent (length, count)\tcategory')
    # template for file path
    file_template = '{0}/{1}'

    # save full stats
    stats_outfile = file_template.format(output_directory, 'loci_stats.tsv')
    lines = [header]
    for locus, values in loci_stats.items():
        most_frequent_str = ';'.join([str(v) for v in values[4]])
        str_values = [values[5]] + [str(v) for v in values[0:4]] \
            + [most_frequent_str] + [values[7]]
        lines.append('\t'.join(str_values))

    text = '\n'.join(lines)
    write_text(text, stats_outfile)

    # create HTML with subplots for loci with alleles outside length threshold
    # and loci with multiple mode values
    # get list of loci
    subplot_traces = [(v[5], v[6]) for k, v in loci_stats.items()
                      if v[6] is not None]

    if len(subplot_traces) > 0:
        # create Figure object for subplots
        subplot_fig = create_subplots(len(subplot_traces),
                                      2,
                                      [v[0] for v in subplot_traces])
        subplot_fig = add_multiple_traces([v[1] for v in subplot_traces],
                                          subplot_fig,
                                          2)

        # adjust figure height and layout
        fig_height = (statistics.math.ceil(len(subplot_fig.data)/2))*300
        subplot_fig.update_layout(title='Length distribution for loci that '
                                        'contain alleles 20% smaller or '
                                        'larger than the allele length mode '
                                        'and/or multimodal distributions.',
                                  height=fig_height,
                                  showlegend=False,
                                  bargap=0.1,
                                  template='ggplot2')

        # log transform axes for more compact scaling
        subplot_fig.update_xaxes(type='log',
                                 title_text='Allele length (log)')
        subplot_fig.update_yaxes(type='log',
                                 title_text='Number of alleles (log)')

        subplot_fig_plotfile = file_template.format(output_directory,
                                                    'loci_barplots.html')
        plot(subplot_fig,
             filename=subplot_fig_plotfile,
             auto_open=False)


def parse_arguments():

    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument('-i', '--schema-directory', type=str,
                        required=True,
                        dest='schema_directory',
                        help='Path to the schema\'s directory.')

    parser.add_argument('-o', '--output-directory', type=str,
                        required=True,
                        dest='output_directory',
                        help='Path to the output directory.')

    parser.add_argument('-s', '--size-threshold', type=float,
                        required=False, default=0.2,
                        dest='size_threshold',
                        help='Locus length mode variation '
                             'threshold for outliers.')

    parser.add_argument('-ft', '--frequency-threshold', type=float,
                        required=False, default=0.3,
                        dest='frequency_threshold',
                        help='')

    args = parser.parse_args()

    return args


if __name__ == '__main__':

    args = parse_arguments()
    main(**vars(args))
