import os
import pandas as pd
import plotly
import plotly.graph_objects as go
import numpy as np

def render_line_plot(file_path, output_path, columns, ascending):
    """
    Render line plot for each column name in input.
    
    Parameters
    ----------
    file_path : str
        Path to the file to read.
    output_path : str
        Where to write HTML file.
    columns : list
        List that contains the columns ids to create line plot.
    ascending : bool
        Bool to choose order to sort (False for highest to smallest).
    
    Returns
    -------
    Writes HTML files inside the output path.
    """
    
    df = pd.read_csv(file_path, sep="\t")
    for column_id in columns:
        fig = go.Figure(data=go.Scatter(y=df[column_id].sort_values(ascending=ascending)))
        html_path = os.path.join(output_path, f'line_plot_{column_id}.html')
        fig.write_html(html_path)
        
def render_histogram(file_path, output_path, columns):
    """
    Render line plot for each column name in input.
    
    Parameters
    ----------
    file_path : str
        Path to the file to read.
    output_path : str
        Where to write HTML file.
    columns : list
        List that contains the columns ids to create histogram plot.
    
    Returns
    -------
    Writes HTML files inside the output path.
    """
    df = pd.read_csv(file_path, sep="\t")
    
    for column_id in columns:
        fig = go.Figure(data=go.Histogram(x=df[column_id]))
        html_path = os.path.join(output_path, f'histogram_{column_id}.html')
        fig.write_html(html_path)
