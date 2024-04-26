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

def create_graph_trace(function, plotname = None, x = None, y = None):
    
    if function == 'boxplot':
        function = go.Box
    elif function == 'scatterplot':
        function = go.Scatter
    elif function == "histogram":
        function = go.Histogram
    
    return function(x= x, y= y, name = plotname)

def generate_plot(traces, title=None, xaxis_title=None, yaxis_title=None):
    """
    Generate a box plot using Plotly.

    Args:
        x_data (list): List of labels for the categories.
        y_data (list of lists): List of lists containing the values for each category.
        title (str, optional): Title of the plot.
        xaxis_title (str, optional): Title of the x-axis.
        yaxis_title (str, optional): Title of the y-axis.

    Returns:
        None
    """

    # Create layout
    layout = go.Layout(title=title, xaxis=dict(title=xaxis_title), yaxis=dict(title=yaxis_title))

    # Create figure
    fig = go.Figure(data=traces, layout=layout)
    
    return fig

def write_fig_to_html(fig, output_path, filename):
    
    html_path = os.path.join(output_path, f'{filename}.html')
    fig.write_html(html_path)
    
def save_plots_to_html(figures, output_path, filename):
    """
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