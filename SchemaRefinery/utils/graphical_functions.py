import os
import plotly.graph_objects as go
import plotly.offline
from plotly.subplots import make_subplots

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

def create_graph_trace(function, x = None , y = None, plotname = None, mode = None):
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
    mode : str, optional
        Mode to create the plot trace (for scatter 'lines' or 'markers')
    
    Returns
    -------
    returns : plotly.graph_objects object
        Returns a trace.

    """
    if function == 'scatter':
        function = go.Scatter
    elif function == 'boxplot':
        function = go.Box
    elif function == 'scatterplot':
        function = go.Scatter
    elif function == 'histogram':
        function = go.Histogram
    elif function == 'violin':
        function = go.Violin
    
    return function(x= x, y= y, name = plotname, mode = mode)

def generate_plot(traces, title=None, xaxis_title=None, yaxis_title=None):
    """
    Generate a plot using Plotly based on input traces.

    Parameters
    ----------
    traces : list
        Contains the list of traces.
    title : str, optional
        Title for the plot.
    xaxis_title : str, optional
        x axis title.
    yaxis_title : str, optional
        y axis title

    Returns
    -------
    fig : plotly.graph_objects object
        Returns a generated plot.
    """

    # Create layout
    layout = go.Layout(title=title,
                       xaxis=dict(title=xaxis_title),
                       yaxis=dict(title=yaxis_title))

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

def plotly_update_layout(fig, title = None, xaxis_title= None, yaxis_title= None, 
                         legend= None, font= None, margin= None, width= None, 
                         height= None, template= None, plot_bgcolor= None, paper_bgcolor= None, 
                         annotations=None, shapes= None, images= None, updatemenus= None, 
                         sliders= None, scene= None, geo= None, showlegend= None):
    """
    Update plotly layout.
    
    Parameters
    ----------

    fig : plotly.graph_objects object
        Plotly go object.
    title : str
    xaxis_title : str
    yaxis_title : str
    legend : dict
        Contains the following keys and values:
            'orientation' : str
                'h' for horizontal, 'v' for vertical.
            'x' : float
                x coordinates of the legend anchor.
            'y' : float
                y coordinates of the legend anchor.
            'xanchor' : str
                Specifies where the legend is anchored.
            'yanchor' : str
                Specifies where the legend is anchored.
    font : dict
        Contains the following keys and values:
            'family' : str
                Font family.
            'size' : int
                Font size.
            'color' : str
                Font color.
    margin : dict
        Contains the following keys and values:
            'l' : int
                Left margin.
            'r' : int
                Right margin.
            't' : int
                Top margin.
            'b' : int
                Bottom margin.
    width : int
        Width of the plot.
    height : int
        Height of the plot.
    template : str or dict
        Plotly template name or template object to apply to the plot.
    plot_bgcolor : str
        Background color of the plot.
    paper_bgcolor : str
        Background color of the plot paper (the area outside the plot).
    annotations : list
        List of dictionaries specifying annotations to be added to the plot.
    shapes : list
        List of dictionaries specifying shapes to be added to the plot.
    images : list
        List of dictionaries specifying images to be added to the plot.
    updatemenus : list
        List of dictionaries specifying update menus to be added to the plot.
    sliders : list
        List of dictionaries specifying sliders to be added to the plot.
    scene : dict
        Dictionary specifying the properties of 3D scenes.
    geo : dict
        Dictionary specifying the properties of geographic maps.
    showlegend : bool
        Whether to show the legend.
        
    Returns
    -------
    fig : plotly.graph_objects object
        Returns a generated plot.
    """
    
    return fig.update_layout(title = title,
                             xaxis_title = xaxis_title,
                             yaxis_title = yaxis_title,
                             legend = legend,
                             font = font,
                             margin = margin,
                             width = width,
                             height = height,
                             template = template,
                             plot_bgcolor = plot_bgcolor,
                             paper_bgcolor = paper_bgcolor,
                             annotations = annotations,
                             shapes = shapes,
                             images = images,
                             updatemenus = updatemenus,
                             sliders = sliders,
                             scene = scene,
                             geo = geo,
                             showlegend = showlegend)
def create_subplots(traces, rows, columns, share_x = None, share_y = None, subplot_titles = None, coords = None):
    """
    Creates subplots from traces using plotly.
    
    Parameters
    ----------
    traces : list
        Contains the list of traces.
    rows : int
        Number of rows.
    columns : int
        Number of columns.
    share_x : bool
        If the traces share same x.
    share_y : bool
        If the traces share same y.
    subplot_titles : list
        Contains the list of subplot titles.
    coords : list
        Contains the list of lists of coords of rows and columns to add each subplot.
    
    fig : plotly.graph_objects object
        Returns subplot.
    """
    fig = make_subplots(rows = rows,
                        cols= columns,
                        shared_xaxes = share_x,
                        shared_yaxes = share_y,
                        subplot_titles = subplot_titles)
    
    for i, trace in enumerate(traces):
        if coords:
            row = traces[i][0]
            column = traces[i][1]
        else:
            row = None
            column = None

        fig.add_trace(trace,
                      row = row,
                      col = column)
    
    return fig

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