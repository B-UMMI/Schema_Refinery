"""
Following script takes and allele call matrix and output several plots and relevant
information about the allele call.

input:
    matrix obtained from allele call.
    
output:
    plots, tables
"""

import os
import argparse
import pandas as pd

import plotly.graph_objs as go
from plotly.offline import plot
from plotly.subplots import make_subplots
import plotly.express as px

import plotly.io as io

io.renderers.default='browser'


def import_matrix(matrix_path):
    
    with open(matrix_path) as profiles:
        
        matrix = pd.read_csv(profiles, sep='\t',low_memory=False)
        
    assembly_ids = matrix.iloc[: , 0]
        
    matrix = matrix.iloc[: , 1:]
    
    classes = ['ALM', 'ASM', 'LNF', 'NIPH',
               'NIPHEM', 'PLOT3', 'PLOT5', 'LOTSC']

    for c in classes:
        if  not matrix[matrix == c].notnull().any().any():
            
            classes.remove(c)

    if matrix.isin(classes).any().any():

        masked_profiles = False
    
    else:
        masked_profiles = True

    return matrix, masked_profiles, assembly_ids, classes

def count_occurences(matrix, masked_profiles, classes):
    
    if not masked_profiles:
        
        return matrix.apply(pd.Series.value_counts).loc[classes]
    
    else:
        pass

def missing_alleles_count_plot(matrix_counts):
    
    matrix_counts_sum = matrix_counts.apply(pd.Series.sum)
    
    matrix_counts_sum = matrix_counts_sum.reset_index()
    
    matrix_counts_sum.columns = ["Locus","Missing"]
    
    histplot_full = px.histogram(matrix_counts_sum,x="Missing",nbins=200,
                        labels={'Missing':'Missing Alleles'}
                        ,text_auto=True)
    
    histplot_partial_100 = px.histogram(matrix_counts_sum[matrix_counts_sum["Missing"] <= 100], x="Missing",nbins=101,
                        labels={'Missing':'Missing Alleles'}
                        ,text_auto=True)
    
    return histplot_full, histplot_partial_100

def stacked_barplot_classes(matrix_counts_t, matrix_counts_sum, assembly_ids):
    
    """
    Takes matrix as input and returns stacked barplot with all classes present
    in the matrix, for each locus.
    
    input:
        matrix with loci as rows and each class as columns, for each loci are 
        counted all assemblies that have been identified as one of each classes.
        
    output:
        stacked barplot for each loci, colored depending on the classses present 
        in the matrix.
    """
    
    matrix_counts_t["Missing_sum"] = matrix_counts_t.sum(axis = 1)
    
    matrix_counts_t = matrix_counts_t.sort_values(by = "Missing_sum")
    
    matrix_counts_t = matrix_counts_t.drop("Missing_sum", axis = 1)
    
    stacked_bar = px.bar(matrix_counts_t)
    
    stacked_bar.update_layout(barmode = "stack",bargap = 0, xaxis_title = "loci", yaxis_title = "Number of assemblies")
    
    stacked_bar.update_xaxes(showticklabels = False)
    
    stacked_bar.add_scatter(x = matrix_counts_sum.sort_values(by="Missing_assemblies").index, 
                            y = matrix_counts_sum.sort_values(by="Missing_assemblies")["Missing_assemblies"], 
                            mode = "lines", line = dict(color="MediumPurple"), name="Sum of all missing alleles")   
    
    return stacked_bar

def main(matrix_path,output_directory):
    
    #import matrix obtained from allele call
    print("loading matrix...")
    matrix , masked_profiles, assembly_ids, classes = import_matrix(matrix_path)
    
    #Count ocurrences of each class in each loci
    print("\ncounting each class occurrence...")
    matrix_counts = count_occurences(matrix, masked_profiles, classes)
    
    #sum  all classes into single category "Missing alleles"
    matrix_counts_sum = matrix_counts.apply(pd.Series.sum).to_frame()
    
    matrix_counts_sum.columns = ["Missing_assemblies"]    
    
    #Transpose matrix
    matrix_counts_t = matrix_counts.transpose()
    
    #Show histogram of the number of loci vs missing alelles
    print("\ngenerating histogram...")
    histplot_full, histplot_partial_100 = missing_alleles_count_plot(matrix_counts)
    
    #stacked barplot of classes for each locus
    print("\ngenerating stacked barplot...")
    stacked_bar = stacked_barplot_classes(matrix_counts_t, matrix_counts_sum, assembly_ids)
    
    #save html
    print("writing html...")
    histplot_full.write_html(os.path.join(output_directory,"histplot_full.html"))

    histplot_partial_100.write_html(os.path.join(output_directory,"histplot_partial_100.html"))
                                    
    stacked_bar.write_html(os.path.join(output_directory,"stacked_barplot.html"))
                           
    #save tables
    print("writing tables...")
    matrix_counts_t.to_csv(os.path.join(output_directory, 'class_table_count.tsv'), sep="\t")
    
                                    
def parse_arguments():

    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument('-i', '--input_matrix', type=str,
                        required=True,
                        dest='matrix_path',
                        help='path to the matrix.')

    parser.add_argument('-o', '--output-directory', type=str,
                        required=True,
                        dest='output_directory',
                        help='Path to the output directory.')


    args = parser.parse_args()

    return args


if __name__ == '__main__':

    args = parse_arguments()
    main(**vars(args))