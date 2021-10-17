#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

"""


import os
import csv
import shutil
import argparse

import numpy as np
import pandas as pd
import mask_matrix as mm
import distance_matrix as dm
from plotly.offline import plot
import Extract_cgAlleles as cg
import plotly.graph_objects as go
import plotly.figure_factory as ff
from scipy.spatial.distance import pdist, squareform
from scipy.cluster.hierarchy import dendrogram, linkage
from matplotlib import pyplot as plt


def simply_return(data):

    return data


def read_tabular(input_file, delimiter='\t'):
    """ Read tabular file.

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


# merge this and next function?
def intra_cluster_stats(dataframe, row_id):
    """
    """

    # get row
    current_row = dataframe.loc[row_id]
    # drop self
    current_row = current_row.drop([row_id])
    # get minimum distance in same cluster
    minimum_distance = current_row.nsmallest(1, keep='all')
    stats_min = (list(minimum_distance.index), minimum_distance[0])
    # get maximum distance in same cluster
    maximum_distance = current_row.nlargest(1, keep='all')
    stats_max = (list(maximum_distance.index), maximum_distance[0])
    # get mean distance in same cluster
    mean_distance = current_row.mean()
    stats_mean = mean_distance

    return [stats_min, stats_max, stats_mean]


def inter_cluster_stats(dataframe, row_id):
    """
    """

    # get minimum distance to strains in other clusters
    minimum_distance = dataframe[row_id].nsmallest(1, keep='all')
    stats_min = (list(minimum_distance.index), minimum_distance[0])
    # get maximum distance to strains in other clusters
    maximum_distance = dataframe[row_id].nlargest(1, keep='all')
    stats_max = (list(maximum_distance.index), maximum_distance[0])
    # get mean distance in same cluster
    mean_distance = dataframe[row_id].mean()
    stats_mean = mean_distance

    return [stats_min, stats_max, stats_mean]


def select_centroids(cluster_stats):
    """
    """

    if len(cluster_stats) > 1:
        means = [(i, j['mean']) for i, j in cluster_stats.items()]
        sorted_means = sorted(means, key=lambda x: x[1])
        centroid = sorted_means[0][0]
    # select singleton as centroid
    else:
        centroid = list(cluster_stats.keys())[0]

    return centroid


# dataframe = distance_matrix
# identifier = k
def centroids_inter_dists(dataframe, centroids_ids, identifier):
    """
    """

    current_cols = ['FILE', identifier]
    df = pd.read_csv(dataframe, sep='\t', usecols=current_cols, index_col='FILE')
    # exclude self and all non-centroids
    outdf = df[df.index.isin(centroids_ids)]
    # get minimum distance to other centroids
    minimum_distance = outdf[identifier].nsmallest(1, keep='all')
    centroid_min = (list(minimum_distance.index), minimum_distance[0])
    # get maximum distance to strains in other clusters
    maximum_distance = outdf[identifier].nlargest(1, keep='all')
    centroid_max = (list(maximum_distance.index), maximum_distance[0])
    # get mean distance in same cluster
    mean_distance = outdf[identifier].mean()
    centroid_mean = mean_distance

    return [centroid_min, centroid_max, centroid_mean]


def compute_stats(matrix, cluster_ids):
    """
    """

    current_cols = ['FILE'] + cluster_ids
    df = pd.read_csv(matrix, sep='\t', usecols=current_cols, index_col='FILE')

    cluster_stats = {e: {} for e in cluster_ids}

    # do not compute statistics for singletons
    if len(cluster_ids) > 1:
        # compute statistics for the cluster
        for e in cluster_ids:
            strain_stats = intra_cluster_stats(df, e)
            cluster_stats[e]['min'] = strain_stats[0]
            cluster_stats[e]['max'] = strain_stats[1]
            cluster_stats[e]['mean'] = strain_stats[2]
    # add empty string to singletons
    else:
        for e in cluster_ids:
            cluster_stats[e]['min'] = ([''], '')
            cluster_stats[e]['max'] = ([''], '')
            cluster_stats[e]['mean'] = ''

    # determine centroid candidates based on minimum mean intra-cluster distance
    cluster_centroid = select_centroids(cluster_stats)
    cluster_stats['centroids'] = cluster_centroid

    return cluster_stats


# function to cluster based on distance matrix with selected algorithm from scipy
def clusters_stats(distance_matrix, clusters):
    """
    """

    # use Pandas to get columns for each cluster
    stats = {}
    for k, v in clusters.items():
        cluster_stats = compute_stats(distance_matrix, v)
        stats[k] = cluster_stats

    return stats


def cluster_strains(distance_matrix, mode):
    """
    """

    df = pd.read_csv(distance_matrix, sep='\t', index_col='FILE')
    df_upper = np.triu(df, k=0)
    df_upper = pd.DataFrame(df_upper, columns=df.columns, index=df.index)

    clusters = linkage(df_upper)


def clusters_boxplot(distance_matrix, clusters):
    """
    """

    colors = ['#f7fcf0', '#ccebc5', '#4eb3d3',
              '#0868ac']

    # use Pandas to get columns for each cluster
    i = 0
    traces = []
    for k, v in clusters.items():
        current_cols = ['FILE'] + v
        df = pd.read_csv(distance_matrix, sep='\t',
                         usecols=current_cols, index_col='FILE')

        # get lines with cluster samples
        df = df.loc[current_cols[1:]]

        # get distance values
        cluster_array = df.to_numpy()
        upper_values = cluster_array[np.triu_indices(len(v), k=1)]
        values = upper_values.tolist()

        # create trace for cluster's boxplot
        trace = go.Box(y=values, name=k,
                       marker=dict(color=colors[i],
                                   line=dict(width=1, color='#252525')),
                       fillcolor=colors[i],
                       line_color='#252525',
                       boxpoints='outliers', jitter=0.5)

        traces.append(trace)

        if i == (len(colors)-1):
            i = 0
        else:
            i += 1

    return traces


def clusters_heatmap(distance_matrix):
    """
    """

    df = pd.read_csv(distance_matrix, sep='\t', index_col='FILE')

    # get distance values
    cluster_array = df.to_numpy()
    values = cluster_array.tolist()
    # reorder matrix lines to group ids in same clusters!

    matrix_header = list(df.columns)

    heatmap = go.Heatmap(z=values[::-1],
                         x=matrix_header,
                         y=matrix_header[::-1],
                         colorscale='Viridis')

    return heatmap


def cluster_dendogram(distance_matrix_file, output_file,
                          linkage_function, distance_function):
    """
    """

    # Read distance matrix
    distance_df = pd.read_csv(distance_matrix_file,
                              sep='\t',
                              index_col='FILE')

    labels = list(distance_df.index)
    # DF to array
    distance_array = distance_df.to_numpy()
    condensed_array = np.triu(distance_array)

    # Initialize figure by creating upper dendrogram
    fig = ff.create_dendrogram(condensed_array, orientation='bottom',
                               labels=labels, distfun=distance_function,
                               linkagefun=linkage_function)

    fig.for_each_trace(lambda trace: trace.update(visible=False))

    # change yaxis for dendogram traces
    for i in range(len(fig['data'])):
        fig['data'][i]['yaxis'] = 'y2'

    # Create Side Dendrogram
    dendro_side = ff.create_dendrogram(condensed_array, orientation='right',
                                       distfun=distance_function,
                                       linkagefun=linkage_function)

    # add traces from right dendogram to existing figure object
    for i in range(len(dendro_side['data'])):
        dendro_side['data'][i]['xaxis'] = 'x2'

    # Add Side Dendrogram Data to Figure
    for data in dendro_side['data']:
        fig.add_trace(data)

    # Create Heatmap
    dendro_leaves = dendro_side['layout']['yaxis']['ticktext']
    dendro_leaves = list(map(int, dendro_leaves))
    heat_data = distance_array
    # reorder data points so that they align with the  dendogram leaves
    # align with right dendogram
    heat_data = heat_data[dendro_leaves,:]
    # align with bottom dendogram
    heat_data = heat_data[:,dendro_leaves]

    heatmap = [go.Heatmap(x = dendro_leaves,
                          y = dendro_leaves,
                          z = heat_data,
                          colorscale = 'Blues')]

    heatmap[0]['x'] = fig['layout']['xaxis']['tickvals']
    heatmap[0]['y'] = dendro_side['layout']['yaxis']['tickvals']

    # Add Heatmap Data to Figure
    for data in heatmap:
        fig.add_trace(data)

    # Edit Layout
    fig.update_layout({'width': 1200, 'height': 900,
                       'showlegend':False, 'hovermode': 'closest'})

    # Edit xaxis with heatmap
    fig.update_layout(xaxis={'domain': [.15, 1],
                                      'mirror': False,
                                      'showgrid': False,
                                      'showline': False,
                                      'zeroline': False,
                                      'showticklabels': False,
                                      'ticks':''})

    # Edit xaxis2 with dendogram
    fig.update_layout(xaxis2={'domain': [0, .15],
                                       'mirror': False,
                                       'showgrid': False,
                                       'showline': False,
                                       'zeroline': False,
                                       'showticklabels': False,
                                       'ticks':''})

    # Edit yaxis with heatmap and side dendogram
    fig.update_layout(yaxis={'domain': [0, 1],
                                      'mirror': False,
                                      'showgrid': False,
                                      'showline': False,
                                      'zeroline': False,
                                      'showticklabels': False,
                                      'ticks': ''})

    # Edit yaxis2 with top dendogram
    fig.update_layout(yaxis2={'domain':[.825, .975],
                                       'mirror': False,
                                       'showgrid': False,
                                       'showline': False,
                                       'zeroline': False,
                                       'showticklabels': False,
                                       'ticks':''})

    fig.update_layout(template='seaborn', paper_bgcolor='rgba(0,0,0,0)',
                      plot_bgcolor='rgba(0,0,0,0)',
                      title='Dendogram and allelic differences heatmap')

    plot(fig, filename=output_file, auto_open=False)


allelecall_results = '/home/rfm/Desktop/rfm/Lab_Analyses/GAS_PrepExternalSchema/cluster_distances/Coelho/coelho_allelecall_results_clean.tsv'
output_directory = '/home/rfm/Desktop/rfm/Lab_Analyses/GAS_PrepExternalSchema/cluster_distances/Coelho/test_results'
clusters = '/home/rfm/Desktop/rfm/Lab_Analyses/GAS_PrepExternalSchema/cluster_distances/Coelho/coelho_clusters.tsv'
cpu_cores = 1
def main(allelecall_results, output_directory, clusters, cpu_cores):

    # Create output directory
    if os.path.isdir(output_directory) is False:
        os.mkdir(output_directory)

    clusters_dict = {}
    if clusters is not None:
        # read clusters    
        clusters_lines = read_tabular(clusters)
        for l in clusters_lines:
            clusters_dict.setdefault(l[1], []).append(l[0])

    # mask matrix
    input_basename = os.path.basename(allelecall_results)
    input_basename = input_basename.split('.tsv')[0]
    masked_results = os.path.join(output_directory, input_basename+'_masked.tsv')
    mm.main(allelecall_results, masked_results, None)

    # import allelecall results
    df = pd.read_csv(masked_results, sep='\t', index_col='FILE')

    clusters_dirs = {}
    first = True
    all_results_dir = os.path.join(output_directory, 'all_results')
    if os.path.isdir(all_results_dir) is False:
        os.mkdir(all_results_dir)
    all_results_file = os.path.join(all_results_dir, 'all_results.tsv')
    clusters_dirs['all_results'] = [all_results_dir, all_results_file]
    for k, v in clusters_dict.items():
        cluster_results = df.loc[v]
        # Create subdir to store cluster results
        cluster_dir = os.path.join(output_directory, k)
        if os.path.isdir(cluster_dir) is False:
            os.mkdir(cluster_dir)

        # save cluster df to cluster dir
        cluster_tsv = os.path.join(cluster_dir, '{0}_allelecall_results.tsv'.format(k))
        cluster_results.to_csv(cluster_tsv, sep='\t')

        clusters_dirs[k] = [cluster_dir, cluster_tsv]

        # append cluster results to file with all results
        if first is True:
            cluster_results.to_csv(all_results_file, sep='\t', mode='a')
            first = False
        else:
            cluster_results.to_csv(all_results_file, sep='\t', mode='a', header=False)

    # determine cgMLST for all clusters
    for k, v in clusters_dirs.items():
        cg_dir = os.path.join(v[0], '{0}_cgMLST'.format(k))
        current_cg = cg.main(v[1], cg_dir, 1.0, [], [])
        cg_file = os.path.join(cg_dir, 'cgMLST.tsv')
        clusters_dirs[k].append(cg_file)

    # create directories to store files with distances
    for k, data in clusters_dirs.items():
        # wgMLST
        wgMLST_distance_dir = os.path.join(data[0], 'wgMLST_{0}_distances'.format(k))
        os.mkdir(wgMLST_distance_dir)

        # cgMLST
        cgMLST_distance_dir = os.path.join(data[0], 'cgMLST_{0}_distances'.format(k))
        os.mkdir(cgMLST_distance_dir)

        clusters_dirs[k].extend([wgMLST_distance_dir, cgMLST_distance_dir])

    # determine distance matrices
    processed = 0
    ad_template = '{0}_allelic_differences.tsv'
    sl_template = '{0}_shared_loci.tsv'
    for k, v in clusters_dirs.items():
        input_basename = os.path.basename(v[1])
        # remove extension that is after last '.'
        input_basename = '.'.join(input_basename.split('.')[0:-1])

        # wgMLST
        tmp_directory = os.path.join(v[3], 'tmp')
        os.mkdir(tmp_directory)

        genome_ids = dm.get_sample_ids(v[1], delimiter='\t')

        np_matrix = dm.tsv_to_nparray(v[1])
        rows_indexes = [i for i in range(len(np_matrix))]
        results = dm.compute_distances(rows_indexes, np_matrix, genome_ids, tmp_directory)

        #merged = dm.merge_dictionaries([results])

        print('\nCreating distance matrix...', end='')
        # create files with headers
        col_ids = ['FILE'] + genome_ids
        output_pairwise = os.path.join(v[3],
                                       ad_template.format(input_basename))
        output_p = os.path.join(v[3],
                                sl_template.format(input_basename))

        # import arrays per genome and save to matrix file
        results2 = dm.write_matrices(results, genome_ids,
                                     output_pairwise, output_p, col_ids)

        # add 1 to include header
        symmetric_allelic_differences = dm.symmetrify_matrix(output_pairwise,
                                                             len(genome_ids)+1,
                                                             tmp_directory)
        symmetric_shared_loci = dm.symmetrify_matrix(output_p,
                                                     len(genome_ids)+1,
                                                     tmp_directory)

        shutil.rmtree(tmp_directory)
        clusters_dirs[k].append(symmetric_allelic_differences)

        # cgMLST
        input_basename = os.path.basename(v[2])
        # remove extension that is after last '.'
        input_basename = '.'.join(input_basename.split('.')[0:-1])
        tmp_directory = os.path.join(v[4], 'tmp')
        os.mkdir(tmp_directory)

        np_matrix = dm.tsv_to_nparray(v[2])
        rows_indexes = [i for i in range(len(np_matrix))]
        # running like this seems to work without ever increasing memory usage...
        results = dm.compute_distances(rows_indexes, np_matrix,
                                       genome_ids, tmp_directory)

        print('\nCreating distance matrix...', end='')
        # create files with headers
        output_pairwise = os.path.join(v[4],
                                       ad_template.format(input_basename))
        output_p = os.path.join(v[4],
                                sl_template.format(input_basename))

        # import arrays per genome and save to matrix file
        results2 = dm.write_matrices(results, genome_ids,
                                    output_pairwise, output_p, col_ids)

        # add 1 to include header
        symmetric_allelic_differences = dm.symmetrify_matrix(output_pairwise,
                                                             len(genome_ids)+1,
                                                             tmp_directory)
        symmetric_shared_loci = dm.symmetrify_matrix(output_p,
                                                     len(genome_ids)+1,
                                                     tmp_directory)

        processed += 1
        print('Processed {0}'.format(processed))

        shutil.rmtree(tmp_directory)
        clusters_dirs[k].append(symmetric_allelic_differences)

###################################

    # running the main function from the distance_matrix.py script originates a memory leak
    # could not determine the cause
    # running other functions that are not the main function does not cause a memory leak

###################################

    # determine intra-cluster stats for all clusters
    all_stats = {}
    for k, v in clusters_dirs.items():
        # determine wg and cg stats for all samples
        if k == 'all_results':
            current_clusters = clusters_dict
        # determine wg and cg stats per cluster
        else:
            current_clusters = {k: clusters_dict[k]}

        # wgMLST
        wgMLST_file = v[5]
        wg_stats = clusters_stats(wgMLST_file, current_clusters)

        # cgMLST
        cgMLST_file = v[6]
        cg_stats = clusters_stats(cgMLST_file, current_clusters)

        all_stats[k] = [wg_stats, cg_stats]

    # create plots for all clusters
    traces = {}
    for k, v in clusters_dirs.items():
        if k == 'all_results':
            current_clusters = clusters_dict
        else:
            current_clusters = {k: clusters_dict[k]}

        # wgMLST
        wgMLST_file = v[5]
        wg_boxplot_trace = clusters_boxplot(wgMLST_file, current_clusters)
        wg_heatmap_trace = clusters_heatmap(wgMLST_file)

        # cgMLST
        cgMLST_file = v[6]
        cg_boxplot_trace = clusters_boxplot(cgMLST_file, current_clusters)
        cg_heatmap_trace = clusters_heatmap(cgMLST_file)

        traces[k] = [[wg_boxplot_trace, wg_heatmap_trace],
                     [cg_boxplot_trace, cg_heatmap_trace]]

    # determine distances between cluster centroids
    # get centroids
    # wgMLST
    wgMLST_centroids = {v['centroids']: k
                        for k, v in all_stats['all_results'][0].items()}

    # only compute if there is more than one cluster
    if len(wgMLST_centroids) > 1:
        distance_matrix = clusters_dirs['all_results'][5]
        for k, v in wgMLST_centroids.items():
            centroids_ids = list(wgMLST_centroids.keys())
            centroids_ids.remove(k)
            centroid_stats = centroids_inter_dists(distance_matrix, centroids_ids, k)
            all_stats['all_results'][0][v]['centroid_min'] = centroid_stats[0]
            all_stats['all_results'][0][v]['centroid_max'] = centroid_stats[1]
            all_stats['all_results'][0][v]['centroid_mean'] = centroid_stats[2]

    # cgMLST
    cgMLST_centroids = {v['centroids']: k
                        for k, v in all_stats['all_results'][1].items()}

    # only compute if there is more than one cluster
    if len(cgMLST_centroids) > 1:
        distance_matrix = clusters_dirs['all_results'][6]
        for k, v in cgMLST_centroids.items():
            centroids_ids = list(cgMLST_centroids.keys())
            centroids_ids.remove(k)
            centroid_stats = centroids_inter_dists(distance_matrix, centroids_ids, k)
            all_stats['all_results'][1][v]['centroid_min'] = centroid_stats[0]
            all_stats['all_results'][1][v]['centroid_max'] = centroid_stats[1]
            all_stats['all_results'][1][v]['centroid_mean'] = centroid_stats[2]

    # get inter-cluster values for all samples
    wgMLST_file = clusters_dirs['all_results'][5]
    cgMLST_file = clusters_dirs['all_results'][6]
    clusters = [k for k in all_stats if k != 'all_results']
    for c in clusters:
        # wgMLST
        current_data = all_stats['all_results'][0][c]
        cluster_ids = [k for k in current_data if 'centroid' not in k]
        current_cols = ['FILE'] + cluster_ids

        wgMLST_df = pd.read_csv(clusters_dirs['all_results'][5],
                                usecols=current_cols,
                                sep='\t',
                                index_col='FILE')

        outdf = wgMLST_df[~wgMLST_df.index.isin(cluster_ids)]
        # do not compute values if there is only one cluster
        if len(outdf) != 0:
            for e in cluster_ids:
                strain_stats = inter_cluster_stats(outdf, e)
                all_stats['all_results'][0][c][e]['min_out'] = strain_stats[0]
                all_stats['all_results'][0][c][e]['max_out'] = strain_stats[1]
                all_stats['all_results'][0][c][e]['mean_out'] = strain_stats[2]

        # cgMLST
        cgMLST_df = pd.read_csv(clusters_dirs['all_results'][6],
                                usecols=current_cols,
                                sep='\t',
                                index_col='FILE')

        outdf = cgMLST_df[~cgMLST_df.index.isin(cluster_ids)]
        # do not compute values if there is only one cluster
        if len(outdf) != 0:
            for e in cluster_ids:
                strain_stats = inter_cluster_stats(outdf, e)
                all_stats['all_results'][1][c][e]['min_out'] = strain_stats[0]
                all_stats['all_results'][1][c][e]['max_out'] = strain_stats[1]
                all_stats['all_results'][1][c][e]['mean_out'] = strain_stats[2]

    # Create HTML files with plots
    for k, v in traces.items():
        outdir = os.path.join(output_directory, k)

        # Create HTML with Boxplots
        # wgMLST
        data = v
        fig = go.Figure()
        wgMLST_html = os.path.join(outdir, 'wgMLST_clusters_boxplots.html')
        for t in data[0][0]:
            fig.add_trace(t)

        fig.update_layout(template='seaborn', title='Intra-emm comparison (number of allelic differences)',
                          xaxis_title='',
                          yaxis_title='Allelic differences',
                          legend_title='Clusters',
                          font=dict(size=14))

        plot(fig, filename=wgMLST_html, auto_open=False)

        # cgMLST
        fig = go.Figure()
        cgMLST_html = os.path.join(outdir, 'cgMLST_clusters_boxplots.html')
        for t in data[1][0]:
            fig.add_trace(t)

        fig.update_layout(template='seaborn', title='Intra-emm comparison (number of allelic differences)',
                          xaxis_title='',
                          yaxis_title='Allelic differences',
                          legend_title='Clusters',
                          font=dict(size=14))

        plot(fig, filename=cgMLST_html, auto_open=False)

        # Heatmap
        # wgMLST
        fig = go.Figure()
        wgMLST_html = os.path.join(outdir, 'wgMLST_heatmap.html')
        fig.add_trace(data[0][1])
        fig.update_layout(template='seaborn', title='Distance matrix (number of allelic differences)')

        plot(fig, filename=wgMLST_html, auto_open=False)

        # cgMLST
        fig = go.Figure()
        cgMLST_html = os.path.join(outdir, 'cgMLST_heatmap.html')
        fig.add_trace(data[1][1])
        fig.update_layout(template='seaborn', title='Distance matrix (number of allelic differences)')

        plot(fig, filename=cgMLST_html, auto_open=False)

    # Create plot with boxplots for the separate analysis of each cluster
    # this is different than the boxplots for the analysis of the complete results
    # these boxplots cannot be directly compared because the set of loci in each might be different
    wgMLST_fig = go.Figure()
    wgMLST_html = os.path.join(output_directory, 'wgMLST_clusters_boxplots.html')
    cgMLST_fig = go.Figure()
    cgMLST_html = os.path.join(output_directory, 'cgMLST_clusters_boxplots.html')
    for k, v in traces.items():
        if k != 'all_results':
            wgMLST_fig.add_trace(v[0][0][0])
            cgMLST_fig.add_trace(v[1][0][0])

    plot(wgMLST_fig, filename=wgMLST_html, auto_open=False)
    plot(cgMLST_fig, filename=cgMLST_html, auto_open=False)

    # Create output files with stats
    file_headers = ['FILE', 'min_distance', 'min_ID', 'max_distance',
                    'max_ID', 'mean_distance', 'min_out_distance',
                    'min_out_ID', 'max_out_distance', 'max_out_ID',
                    'mean_out_distance', 'cluster']
    cluster_headers = file_headers[0:6]
    for k, v in all_stats.items():
        if k != 'all_results':
            # wgMLST
            wgMLST_stats = v[0][k]
            outfile = os.path.join(output_directory, k, 'wgMLST_stats.tsv')
            wgMLST_lines = [cluster_headers]
            for h, j in wgMLST_stats.items():
                if 'centroid' not in h:
                    current_id = h
                    min_distance = j['min']
                    min_id = min_distance[0][0]
                    min_distance = str(min_distance[1])

                    max_distance = j['max']
                    max_id = max_distance[0][0]
                    max_distance = str(max_distance[1])

                    mean_distance = j['mean']
                    if mean_distance != '':
                        mean_distance = str(round(mean_distance, 2))

                    wgMLST_lines.append([h, min_distance, min_id,
                                         max_distance, max_id,
                                         mean_distance])

            wgMLST_lines = ['\t'.join(l) for l in wgMLST_lines]
            wgMLST_text = '\n'.join(wgMLST_lines)
            with open(outfile, 'w') as outh:
                outh.write(wgMLST_text+'\n')

            # cgMLST
            cgMLST_stats = v[1][k]
            outfile = os.path.join(output_directory, k, 'cgMLST_stats.tsv')
            cgMLST_lines = [cluster_headers]
            for h, j in cgMLST_stats.items():
                if 'centroid' not in h:
                    current_id = h
                    min_distance = j['min']
                    min_id = min_distance[0][0]
                    min_distance = str(min_distance[1])

                    max_distance = j['max']
                    max_id = max_distance[0][0]
                    max_distance = str(max_distance[1])

                    mean_distance = j['mean']
                    if mean_distance != '':
                        mean_distance = str(round(mean_distance, 2))

                    cgMLST_lines.append([h, min_distance, min_id,
                                         max_distance, max_id,
                                         mean_distance])

            cgMLST_lines = ['\t'.join(l) for l in cgMLST_lines]
            cgMLST_text = '\n'.join(cgMLST_lines)
            with open(outfile, 'w') as outh:
                outh.write(cgMLST_text+'\n')

    # Create output files with global stats
    global_stats = all_stats['all_results']
    clusters_ids = list(global_stats[0].keys())
    wgMLST_lines = [file_headers]
    cgMLST_lines = [file_headers]
    for k in clusters_ids:
        wgMLST_stats = global_stats[0][k]
        for h, j in wgMLST_stats.items():
            if 'centroid' not in h:
                current_id = h
                min_distance = j['min']
                min_id = min_distance[0][0]
                min_distance = str(min_distance[1])

                max_distance = j['max']
                max_id = max_distance[0][0]
                max_distance = str(max_distance[1])

                mean_distance = j['mean']
                if mean_distance != '':
                    mean_distance = str(round(mean_distance, 2))

                min_distance_out = j['min_out']
                min_id_out = min_distance_out[0][0]
                min_distance_out = str(min_distance_out[1])

                max_distance_out = j['max_out']
                max_id_out = min_distance_out[0][0]
                max_distance_out = str(max_distance_out[1])

                mean_distance_out = str(round(j['mean_out'], 2))

                wgMLST_lines.append([h, min_distance, min_id,
                                     max_distance, max_id,
                                     mean_distance, min_distance_out,
                                     min_id_out, max_distance_out,
                                     max_id_out, mean_distance_out,
                                     k])

        cgMLST_stats = global_stats[1][k]
        outfile = os.path.join(output_directory, 'all_results', 'cgMLST_stats.tsv')
        for h, j in cgMLST_stats.items():
            if 'centroid' not in h:
                current_id = h
                min_distance = j['min']
                min_id = min_distance[0][0]
                min_distance = str(min_distance[1])

                max_distance = j['max']
                max_id = max_distance[0][0]
                max_distance = str(max_distance[1])

                mean_distance = j['mean']
                if mean_distance != '':
                    mean_distance = str(round(mean_distance, 2))

                min_distance_out = j['min_out']
                min_id_out = min_distance_out[0][0]
                min_distance_out = str(min_distance_out[1])

                max_distance_out = j['max_out']
                max_id_out = min_distance_out[0][0]
                max_distance_out = str(max_distance_out[1])

                mean_distance_out = str(round(j['mean_out'], 2))

                cgMLST_lines.append([h, min_distance, min_id,
                                     max_distance, max_id,
                                     mean_distance, min_distance_out,
                                     min_id_out, max_distance_out,
                                     max_id_out, mean_distance_out,
                                     k])

    outfile = os.path.join(output_directory, 'all_results', 'wgMLST_stats.tsv')
    with open(outfile, 'w') as outh:
        wgMLST_lines = ['\t'.join(l) for l in wgMLST_lines]
        wgMLST_text = '\n'.join(wgMLST_lines)
        outh.write(wgMLST_text+'\n')

    outfile = os.path.join(output_directory, 'all_results', 'cgMLST_stats.tsv')
    with open(outfile, 'w') as outh:
        cgMLST_lines = ['\t'.join(l) for l in cgMLST_lines]
        cgMLST_text = '\n'.join(cgMLST_lines)
        outh.write(cgMLST_text+'\n')

    #######################
    # Dendogram and heatmap

    # define linkage function to use
    # Perform hierarchical/agglomerative clustering
    # single/min/nearest linkage on the condensed distance matrix
    linkage_function = lambda x: linkage(x, 'single')

    # distfun simply reads the distance matrix that was computed
    # the implementation of ff.create_dendrogram applied the distfun
    # to the whole array. By defining a distfun that returns the passed
    # argument we can pass a precomputed distance matrix
    distance_function = simply_return

    for k, v in clusters_dirs.items():
        # wgMLST
        distance_matrix_file = v[5]
        output_file = os.path.join(v[0], 'wgMLST_dendogram.html')
        cluster_dendogram(distance_matrix_file, output_file,
                          linkage_function, distance_function)

        # cgMLST
        distance_matrix_file = v[6]
        output_file = os.path.join(v[0], 'cgMLST_dendogram.html')
        cluster_dendogram(distance_matrix_file, output_file,
                          linkage_function, distance_function)


def parse_arguments():

    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument('-d', '--allelecall-results', type=str,
                        required=True, dest='allelecall_results',
                        help='')
    
    parser.add_argument('-o', '--output-directory', type=str,
                        required=True, dest='output_directory',
                        help='')

    parser.add_argument('-m', '--mode', type=str,
                        required=True, dest='mode',
                        help='')

    parser.add_argument('-o', '--output-directory', type=str,
                        required=True, dest='output_directory',
                        help='')

    parser.add_argument('--clusters', type=str,
                        required=False, dest='clusters',
                        help='')

    parser.add_argument('-c', '--cpu-cores', type=int,
                        required=False, default=1,
                        dest='cpu_cores',
                        help='')

    args = parser.parse_args()

    return args


if __name__ == "__main__":

    args = parse_arguments()
    main(**vars(args))
