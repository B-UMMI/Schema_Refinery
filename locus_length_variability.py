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
import shutil
import argparse
import csv
import statistics
import hashlib
from collections import Counter

from Bio import SeqIO
import plotly.graph_objs as go
from plotly.offline import plot
from plotly.subplots import make_subplots

RED = '#bd0026'
YELLOW = '#fec44f'
BLUE = '#0868ac'
GREEN = '#238b45'
GRAY = '#969696'

def checkAndMakeDirectory(outdir: str):
    if not os.path.isdir(outdir):
        os.mkdir(outdir)

def listTSVFile(file: str):
    with open(file, 'r') as table1:
        lines_table_1 = list(csv.reader(table1, delimiter='\t'))
    return lines_table_1

def exportDictToTSVFile(filename: str, dict: dict, outdir: str, first_row: list):

    csv_file = f"{filename}.csv"

    rows = []

    for hash, ids_list in dict.items():
        rows.append([hash, " | ".join(ids_list)])

    with open(os.path.join(outdir, csv_file), 'w') as csvfile: 
        # creating a csv writer object 
        csvwriter = csv.writer(csvfile, delimiter='\t') 

        # writing the data rows 
        csvwriter.writerow(first_row) #put the collumn names
        csvwriter.writerows(rows)

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


def histogram_trace(xdata, ydata, marker_colors, name):
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
                   marker_color=marker_colors,
                   name=name,
                   hoverlabel=dict(namelength = -1))

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
        for bar in t:
            fig.append_trace(bar, row_start, col_start)
        if col_start < max_cols:
            col_start += 1
        else:
            col_start = 1
            row_start += 1

    return fig

def add_multiple_traces_paralogous(traces_dict: dict, fig, max_cols, colors):
    row_start = 1
    col_start = 1
    for t in traces_dict:
        for bar in t[1]:
            fig.append_trace(bar, row_start, col_start)

        fig.add_vline(x=colors[t[0]][GREEN][0], row=row_start, col=col_start, line_width=3, line_dash="dash", line_color='black', opacity=1)
        if colors[t[0]][BLUE]:
            for x_val in colors[t[0]][BLUE]:
                fig.add_vline(x=x_val, row=row_start, col=col_start, line_width=1, line_color='blue', opacity=1)
        if colors[t[0]][RED]:
            fig.add_vrect(x0=min(colors[t[0]][RED]), x1=max(colors[t[0]][RED]), row=row_start, col=col_start, fillcolor="red", opacity=0.2, layer="below", line_width=0)
        if colors[t[0]][YELLOW]:
            fig.add_vrect(x0=min(colors[t[0]][YELLOW]), x1=max(colors[t[0]][YELLOW]), row=row_start, col=col_start, fillcolor="yellow", opacity=0.5, layer="below", line_width=0)

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

def parse_paralogous(join_paralogous_file: str, schema_directory: str, output_directory: str):
    all_lines = listTSVFile(join_paralogous_file)
    all_lines.pop(0)

    hash_paralogous_correspondence = {}
    all_files = set()

    checkAndMakeDirectory(os.path.join(schema_directory, "joined_paralogous"))

    for line in all_lines:
        ids_list = line[0].split(' | ')
        ids_paths = [os.path.join(schema_directory, f"{id}.fasta") for id in ids_list]

        hash_id = hashlib.md5(line[0].encode()).hexdigest()

        hash_filename = hashlib.md5(line[0].encode()).hexdigest()

        new_file_path = os.path.join(schema_directory, "joined_paralogous", hash_filename + ".fasta") 
        
        hash_paralogous_correspondence[hash_id] = ids_list

        with open(new_file_path, 'a') as new_file:
            for path in ids_paths:
                with open(path, 'r') as fasta:
                    lines = fasta.readlines()
                    new_file.writelines(lines)

        all_files.add(new_file_path)
        all_files.update(ids_paths)

    exportDictToTSVFile("paralogous_hash_correspondences", hash_paralogous_correspondence, output_directory, ["hash", "file ids"])

    return list(all_files), hash_paralogous_correspondence

def build_locus_graphs(loci_stats: dict, file_template, output_directory: str):
    # create HTML with subplots for loci with alleles outside length threshold
    # and loci with multiple mode values
    # get list of loci
    subplot_traces = [(v[5], v[6]) for k, v in loci_stats.items()
                        if v[6] is not None]

    #print(subplot_traces)
    #print(loci_stats)

    if len(subplot_traces) > 0:
        # create Figure object for subplots
        subplot_fig = create_subplots(len(subplot_traces),
                                      2,
                                      [v[0] for v in subplot_traces])
        subplot_fig = add_multiple_traces([[v[1]] for v in subplot_traces],
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

def build_paralogous_graphs(loci_stats: dict, file_template, output_directory: str, paralogous_groups: dict):
    # create HTML with subplots for loci with alleles outside length threshold
    # and loci with multiple mode values
    # get list of loci

    subplot_traces = [(v[5], v[6]) for k, v in loci_stats.items()]
    colors = {k: v[-1] for k, v in loci_stats.items() if k in paralogous_groups.keys()}

    subplot_traces_hashes = []
    for hash, ids_list in paralogous_groups.items():
        data = []
        for id in ids_list:
            for trace in subplot_traces:
                if trace[0] == id:
                    data.append(trace[1])
        subplot_traces_hashes.append((hash, data))

    #print(colors)

    if len(subplot_traces_hashes) > 0:
        # create Figure object for subplots
        subplot_fig = create_subplots(len(subplot_traces_hashes),
                                      1,
                                      [v[0] for v in subplot_traces_hashes])
        subplot_fig = add_multiple_traces_paralogous([v for v in subplot_traces_hashes],
                                          subplot_fig,
                                          1, colors)

        # adjust figure height and layout
        fig_height = (statistics.math.ceil(len(subplot_fig.data)))*300
        subplot_fig.update_layout(title='Length distribution for loci',
                                  height=fig_height,
                                  showlegend=False,
                                  bargap=0.1,
                                  template='ggplot2',
                                  barmode='stack')

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

def locus_process(schema_directory, output_directory, size_threshold, frequency_threshold):
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
                colors.append(YELLOW)
            elif v >= alm:
                colors.append(RED)
            elif v == most_frequent[0]:
                colors.append(GREEN)
            elif len(high_freq) > 0 and v in [f[0] for f in high_freq]:
                colors.append(BLUE)
            else:
                colors.append(GRAY)

        loci_stats[locus] = [total_alleles, locus_mean,
                             locus_median, locus_mode,
                             frequent_values, locus_id]

        # create histogram for loci that have multimodal distribution or
        # loci that have alleles outside the acceptable length interval
        if locus_mode == 'None' or len(outlier) > 0:
            # create histogram object
            trace = histogram_trace(xdata, ydata, colors, None)
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

    build_locus_graphs(loci_stats, file_template, output_directory)

def paralogous_process(schema_directory, output_directory, size_threshold, frequency_threshold, paralogous_file):
    all_ids_files, paralogous_groups = parse_paralogous(paralogous_file, schema_directory, output_directory)
    schema_loci = [file for file in all_ids_files]

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

        colors = {YELLOW: [], RED: [], GREEN: [], BLUE: []}
        for v in xdata:
            if v <= asm:
                colors[YELLOW].append(v)
            elif v >= alm:
                colors[RED].append(v)
            elif v == most_frequent[0]:
                colors[GREEN].append(v)
            elif len(high_freq) > 0 and v in [f[0] for f in high_freq]:
                colors[BLUE].append(v)

        loci_stats[locus_id] = [total_alleles, locus_mean,
                             locus_median, locus_mode,
                             frequent_values, locus_id]

        # create histogram for loci that have multimodal distribution or
        # loci that have alleles outside the acceptable length interval

        # this trace creation was moved outside of the if to always generate graph information for all loci
        # create histogram object
        trace = histogram_trace(xdata, ydata, None, locus_id)
        if locus_mode == 'None' or len(outlier) > 0:
            if len(outlier) > 0 and locus_mode == 'None':
                loci_stats[locus_id].extend([trace, 'outlier+multimodal'])
            elif len(outlier) > 0:
                loci_stats[locus_id].extend([trace, 'outlier'])
            elif locus_mode == 'None':
                loci_stats[locus_id].extend([trace, 'multimodal'])
        else:
            loci_stats[locus_id].extend([trace, 'single'])

        loci_stats[locus_id].extend([colors])



    # header for file with statistics
    header = ('locus\tnum_alleles\tmean\tmedian\tmode\t'
              'most_frequent (length, count)\tcategory')
    # template for file path
    file_template = '{0}/{1}'


    # build graph html file
    build_paralogous_graphs(loci_stats, file_template, output_directory, paralogous_groups)
    shutil.rmtree(os.path.join(schema_directory, "joined_paralogous"))


    # save full stats to stats file
    stats_outfile = file_template.format(output_directory, 'loci_stats.tsv')
    lines = [header]

    # filter locus from loci_stats
    loci_stats_items = loci_stats.items()
    paralogous_groups_keys = paralogous_groups.keys()
    filtered_loci_stats = {}
    for key, value in loci_stats_items:
        if key in paralogous_groups_keys:
            filtered_loci_stats[key] = value
    
    loci_stats = filtered_loci_stats

    for _, values in loci_stats.items():
        most_frequent_str = ';'.join([str(v) for v in values[4]])
        str_values = [values[5]] + [str(v) for v in values[0:4]] \
            + [most_frequent_str] + [values[7]]
        lines.append('\t'.join(str_values))

    text = '\n'.join(lines)
    write_text(text, stats_outfile)

def main(schema_directory, output_directory, size_threshold, frequency_threshold, paralogous_file):
    checkAndMakeDirectory(output_directory)

    if paralogous_file:
        paralogous_process(schema_directory, output_directory, size_threshold, frequency_threshold, paralogous_file)
    else:
        locus_process(schema_directory, output_directory, size_threshold, frequency_threshold)

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

    parser.add_argument('-pf', '--paralogous-file', type=str,
                        required=False, default='',
                        dest='paralogous_file',
                        help='path to the file with the paralogous ids to be joined')

    args = parser.parse_args()

    return args


if __name__ == '__main__':

    args = parse_arguments()
    main(**vars(args))
