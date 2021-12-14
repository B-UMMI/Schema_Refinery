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
    """ Imports sequences from a FASTA file.

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


def histogram_trace(xdata, ydata, marker_color):
    """ Creates a histogram (go.Bar) trace.

    Parameters
    ----------
    xdata : list
        Sets the x coordinates (sequence
        length values).
    ydata : list
        Sets the y coordinates (counts for
        each sequence length value).
    marker_color : str
        Color to fill the bars.

    Returns
    -------
    trace : plotly.graph_objects.Bar
        Trace with data to display histogram.
    """

    trace = go.Bar(x=xdata,
                   y=ydata,
                   marker_color=marker_color)

    return trace


def create_subplots(nplots, ncols, titles):
    """ Create a Figure objcet with predefined subplots.

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
    """ Add traces to Figure object based on the number
        of cols per row.

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
    """ Writes text to a file.

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


def main(schema_directory, output_directory, size_threshold):

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
        # locus does not have a unique mode
        except Exception as e:
            locus_id = os.path.basename(locus).split('.fasta')[0]
            print(locus_id, '\n', e)
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

        # determine if there are alleles outside (mode +/- (mode*size_threshold))
        asm = most_frequent[0]-(most_frequent[0]*size_threshold)
        alm = most_frequent[0]+(most_frequent[0]*size_threshold)
        outlier = [v
                   for v in freq_ordered[1:]
                   if (v[0] >= alm or v[0] <= asm)]

        xdata = [v[0] for v in freq_ordered]
        ydata = [v[1] for v in freq_ordered]

        loci_stats[locus] = [total_alleles, locus_mean,
                             locus_median, locus_mode,
                             frequent_values, locus_id]

        if locus_mode == 'None' or len(outlier) > 0:
            # create histogram object
            trace = histogram_trace(xdata, ydata, '#0868ac')
            if locus_mode == 'None':
                loci_stats[locus].extend([trace, 'multimodal'])
            elif len(outlier) > 0:
                loci_stats[locus].extend([trace, 'outlier'])
        else:
            loci_stats[locus].extend([None, 'single'])

    # header for file with statistics
    header = 'locus\tnum_alleles\tmean\tmedian\tmode\tmost_frequent\tcategory'
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
    # get list of loci
    outlier_loci = [v[5]
                    for k, v in loci_stats.items()
                    if v[-1] == 'outlier']

    # get traces
    outlier_traces = [v[-2]
                      for k, v in loci_stats.items()
                      if v[-1] == 'outlier']

    # create Figure object for subplots
    outlier_mode_fig = create_subplots(len(outlier_traces), 2, outlier_loci)
    outlier_mode_fig = add_multiple_traces(outlier_traces, outlier_mode_fig,
                                           2)

    # adjust figure height and layout
    fig_height = (statistics.math.ceil(len(outlier_mode_fig.data)/2))*300
    outlier_mode_fig.update_layout(title='Loci with alleles 20% smaller or '
                                         'larger than the allele length mode',
                                   height=fig_height,
                                   showlegend=False,
                                   bargap=0.1,
                                   template='ggplot2')

    # log transform axes for more compact scaling
    outlier_mode_fig.update_xaxes(type='log',
                                  title_text='Allele length (log)')
    outlier_mode_fig.update_yaxes(type='log',
                                  title_text='Number of alleles (log)')

    outlier_mode_plotfile = file_template.format(output_directory,
                                                 'outlier.html')
    plot(outlier_mode_fig, filename=outlier_mode_plotfile, auto_open=False)

    # create HTML with subplots for loci with multiple modes
    multimodal_loci = [v[5]
                       for k, v in loci_stats.items()
                       if v[-1] == 'multimodal']

    multimodal_traces = [v[-2]
                         for k, v in loci_stats.items()
                         if v[-1] == 'multimodal']

    multimodal_fig = create_subplots(len(multimodal_traces), 2, multimodal_loci)
    multimodal_fig = add_multiple_traces(multimodal_traces, multimodal_fig,
                                         2)

    fig_height = (statistics.math.ceil(len(multimodal_fig.data)/2))*300
    multimodal_fig.update_layout(title='Loci with multimodal allele length distributions',
                                 height=fig_height,
                                 showlegend=False,
                                 bargap=0.1,
                                 template='ggplot2')

    multimodal_fig.update_xaxes(type='log', title_text='Allele length (log)')
    multimodal_fig.update_yaxes(type='log', title_text='Number of alleles (log)')

    plot_file = file_template.format(output_directory,
                                     'multimodal.html')
    plot(multimodal_fig, filename=plot_file, auto_open=False)


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

    args = parser.parse_args()

    return args


if __name__ == '__main__':

    args = parse_arguments()
    main(**vars(args))
