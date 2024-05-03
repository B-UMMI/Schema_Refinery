import os
import plotly.graph_objects as go
import plotly.offline

def render_line_chart(df, output_path, columns, labels, ascending):
    """
    Render line chart for each column name in input.
    
    Parameters
    ----------
    df : pandas dataframe
        pandas dataframe
    output_path : str
        Where to write HTML file.
    columns : list
        List that contains the columns ids to create line chart.
    labels : list
        Contains x and y lables to add to the chart.
    ascending : bool
        Bool to choose order to sort (False for highest to smallest).
    
    Returns
    -------
    Writes HTML files inside the output path.
    """
    
    for column_id in columns:
        fig = go.Figure(data=go.Scatter(y=df[column_id].sort_values(ascending=ascending)))
        # Update layout
        fig.update_layout(
        title=f"Line chart of {column_id}",
        xaxis_title="Value",
        yaxis_title="Frequency")
        
        html_path = os.path.join(output_path, f'line_chart_{column_id}.html')
        fig.write_html(html_path)
        
def render_histogram(df, output_path, columns, labels):
    """
    Render line plot for each column name in input.
    
    Parameters
    ----------
    df : pandas dataframe
        pandas dataframe
    output_path : str
        Where to write HTML file.
    columns : list
        List that contains the columns ids to create histogram plot.
    labels : list
        Contains x and y lables to add to the chart.
    
    Returns
    -------
    Writes HTML files inside the output path.
    """
    
    
    for column_id in columns:
        fig = go.Figure(data=go.Histogram(x=df[column_id]))
        # Update layout
        fig.update_layout(
        title=f"Histogram of {column_id}",
        xaxis_title=labels[0],
        yaxis_title=labels[1])
        
        html_path = os.path.join(output_path, f'histogram_{column_id}.html')
        fig.write_html(html_path)

def create_graph_trace(function, x = None , y = None, plotname = None):
    """
    Based on input creates and returns a graph trace using plotly.
    
    Parameters
    ----------
    function : str
        Which plot to trace, can be ['boxplot', 'scatterplot', 'histogram']
    x : pandas dataframe, optional
        Values to add to the trace contained in pandas dataframe.
    y : pandas dataframe, optional
        Values to add to the trace contained in pandas dataframe.
    plotname : str, optional
        The name for the plot.
    
    Returns
    -------
    returns : plotly.graph_objects object
        Returns a trace.

    """
    if function == 'boxplot':
        function = go.Box
    elif function == 'scatterplot':
        function = go.Scatter
    elif function == 'histogram':
        function = go.Histogram
    
    return function(x= x, y= y, name = plotname)

def generate_plot(traces, title=None, xaxis_title=None, yaxis_title=None):
    """
    Generate a plot using Plotly based on input traces.

    Parameters
    ----------
    traces : list
        Contains List of traces.
    title : str, optional
        Title for the plot.
    xaxis_title : str, optional
        x axis title.
    yaxis_title : str, optional
        y axis title

    Returns
    -------
    returns : plotly.graph_objects object
        Returns a generated plot.
    """

    # Create layout
    layout = go.Layout(title=title, xaxis=dict(title=xaxis_title), yaxis=dict(title=yaxis_title))

    # Create figure
    fig = go.Figure(data=traces, layout=layout)
    
    return fig

def write_fig_to_html(fig, output_path, filename):
    """
    Writes the fig plotly object to HTML file.

    Parameters
    ----------
    fig : plotly.graph_objects object
        Figure onject to write to HTML file.
    output_path : str
        Path were to save the HTML file.
    filename : str
        Name for the HTML file.
        
    Returns
    -------
    Writes an HTML file containg the plot to the output directory.
    """
    html_path = os.path.join(output_path, f'{filename}.html')
    fig.write_html(html_path)
    
def save_plots_to_html(figures, output_path, filename):
    """
    Saves a list of plotly.graph_objects object to one HTML file.
    
    Parameters
    ----------
    figures : list
        List that contains the figures to save.
    output_path : str
        Path to the output directory.
    filename : str
        Name for the HTML file.

    Returns
    -------
    Writes an HTML file containg the plots to the output directory.
    """
    
    html_path = os.path.join(output_path, f'{filename}.html')
    # Open a new HTML file
    with open(html_path, 'w') as f:
        # Write HTML header
        f.write('<!DOCTYPE html>\n<html>\n<head><title>Plotly Plots</title></head>\n<body>\n')
    
        # Write each figure's HTML div to the file
        for fig in figures:
            f.write('<div>\n')
            f.write(plotly.offline.plot(fig, output_type='div', include_plotlyjs=True))
            f.write('\n</div>\n')

        # Write HTML footer
        f.write('</body>\n</html>')