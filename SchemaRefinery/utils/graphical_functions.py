import os
import pandas as pd
import plotly.graph_objects as go, Figure
import plotly.offline
from plotly.subplots import make_subplots
from typing import List, Optional, Union, Dict, Any

def render_line_chart(df: pd.DataFrame, output_path: str, columns: List[str], ascending: bool) -> None:
    """
    Render line chart for each column name in input.

    Parameters
    ----------
    df : pd.DataFrame
        Pandas dataframe.
    output_path : str
        Where to write HTML file.
    columns : List[str]
        List that contains the column ids to create line chart.
    ascending : bool
        Bool to choose order to sort (False for highest to smallest).

    Returns
    -------
    None
        Writes HTML files inside the output path.
    """
    for column_id in columns:
        fig: go.Figure = go.Figure(data=go.Scatter(y=df[column_id].sort_values(ascending=ascending)))
        # Update layout
        fig.update_layout(
            title=f"Line chart of {column_id}",
            xaxis_title="Value",
            yaxis_title="Frequency"
        )
        
        html_path = os.path.join(output_path, f'line_chart_{column_id}.html')
        fig.write_html(html_path)

def render_histogram(df: pd.DataFrame, output_path: str, columns: List[str], labels: List[str]) -> None:
    """
    Render histogram for each column name in input.

    Parameters
    ----------
    df : pd.DataFrame
        Pandas dataframe.
    output_path : str
        Where to write HTML file.
    columns : List[str]
        List that contains the column ids to create histogram plot.
    labels : List[str]
        Contains x and y labels to add to the chart.

    Returns
    -------
    None
        Writes HTML files inside the output path.
    """
    for column_id in columns:
        fig: go.Figure = go.Figure(data=go.Histogram(x=df[column_id]))
        # Update layout
        fig.update_layout(
            title=f"Histogram of {column_id}",
            xaxis_title=labels[0],
            yaxis_title=labels[1]
        )
        
        html_path: str = os.path.join(output_path, f'histogram_{column_id}.html')
        fig.write_html(html_path)

def create_graph_trace(
    function: str, 
    x: Optional[Union[pd.DataFrame, List[Any]]] = None, 
    y: Optional[Union[pd.DataFrame, List[Any]]] = None, 
    plotname: Optional[str] = None, 
    mode: Optional[str] = None
) -> go.Figure:
    """
    Based on input creates and returns a graph trace using Plotly.

    Parameters
    ----------
    function : str
        Which plot to trace, can be ['boxplot', 'scatterplot', 'histogram', 'violin'].
    x : Optional[Union[pd.DataFrame, List[Any]]], optional
        Values to add to the trace contained in pandas dataframe or list.
    y : Optional[Union[pd.DataFrame, List[Any]]], optional
        Values to add to the trace contained in pandas dataframe or list.
    plotname : Optional[str], optional
        The name for the plot.
    mode : Optional[str], optional
        Mode to create the plot trace (for scatter 'lines' or 'markers').

    Returns
    -------
    go.Figure
        The Plotly figure object containing the plot.
    """
    if function == 'scatter':
        return go.Scatter(x=x, y=y, name=plotname, mode=mode)
    elif function == 'boxplot':
        function = go.Box
    elif function == 'scatterplot':
        function = go.Scatter
    elif function == 'histogram':
        function = go.Histogram
    elif function == 'violin':
        function = go.Violin
    
    return function(x=x, y=y, name=plotname)

def create_violin_plot(
    y: Optional[Union[List[Any], pd.Series]] = None, 
    x: Optional[Union[str, int, List[Union[str, int]]]] = None, 
    name: Optional[str] = None, 
    orientation: str = 'v', 
    box_visible: bool = True, 
    meanline_visible: bool = False, 
    points: str = 'outliers', 
    line_color: Optional[Union[str, int]] = None, 
    marker: Optional[Dict[str, Any]] = None, 
    opacity: float = 1, 
    side: str = 'both', 
    scalemode: str = 'width', 
    spanmode: str = 'soft', 
    scalegroup: Optional[str] = None
) -> go.Violin:
    """
    Create a violin plot using Plotly.

    Parameters
    ----------
    y : Optional[Union[List[Any], pd.Series]], optional
        The data for the violin plot.
    x : Optional[Union[str, int, List[Union[str, int]]]], optional
        The group labels or positions along the x-axis for the violin plot.
    name : Optional[str], optional
        The given name for the plot.
    orientation : str, optional
        The orientation of the violin plot ('v' for vertical, 'h' for horizontal).
    box_visible : bool, optional
        Whether to display the box plot inside the violin plot (default is True).
    meanline_visible : bool, optional
        Whether to display the mean line inside the violin plot (default is False).
    points : str, optional
        Whether to show points on the violin plot ('outliers', 'suspectedoutliers', or 'all').
    line_color : Optional[Union[str, int]], optional
        The color of the line around the violin plot.
    marker : Optional[Dict[str, Any]], optional
        A dictionary specifying the marker properties for the points (if shown).
    opacity : float, optional
        The opacity of the violin plot (a value between 0 and 1, default is 1).
    side : str, optional
        The side of the violin plot to plot points ('both', 'positive', or 'negative').
    scalemode : str, optional
        The scaling mode for the violin plot ('width' or 'count', default is 'width').
    spanmode : str, optional
        The span mode for the violin plot ('soft' or 'hard', default is 'soft').
    scalegroup : Optional[str], optional
        The name of the group of violins whose widths should be made proportional to the number of samples in each group.

    Returns
    -------
    go.Violin
        The Plotly figure object containing the violin plot.
    """
    return go.Violin(
        y=y, x=x, name=name, orientation=orientation, box_visible=box_visible,
        meanline_visible=meanline_visible, points=points, line_color=line_color,
        marker=marker, opacity=opacity, side=side, scalemode=scalemode,
        spanmode=spanmode, scalegroup=scalegroup
    )

def create_box_plot(
    y: Optional[Union[List[Any], pd.Series]] = None, 
    x: Optional[Union[str, int, List[Union[str, int]]]] = None, 
    name: Optional[str] = None, 
    orientation: str = 'v', 
    boxpoints: str = 'outliers', 
    jitter: float = 0.3
) -> go.Box:
    """
    Create a box plot using Plotly.

    Parameters
    ----------
    y : Optional[Union[List[Any], pd.Series]], optional
        The data for the box plot.
    x : Optional[Union[str, int, List[Union[str, int]]]], optional
        The group labels or positions along the x-axis for the box plot.
    name : Optional[str], optional
        The given name for the plot.
    orientation : str, optional
        The orientation of the box plot ('v' for vertical, 'h' for horizontal).
    boxpoints : str, optional
        Specifies how the data points are displayed ('outliers', 'suspectedoutliers', 'all', or False).
    jitter : float, optional
        Sets the amount of jitter in the box plot points (0 for no jitter, 1 for full jitter).

    Returns
    -------
    go.Box
        The Plotly figure object containing the box plot.
    """
    return go.Box(
        y=y, x=x, name=name, orientation=orientation,
        boxpoints=boxpoints, jitter=jitter
    )

def create_histogram(
    x: Union[List[Any], pd.Series], 
    name: Optional[str] = None, 
    xbins: Optional[Dict[str, Any]] = None, 
    histnorm: Optional[str] = None, 
    orientation: str = 'v'
) -> go.Histogram:
    """
    Create a histogram using Plotly.

    Parameters
    ----------
    x : Union[List[Any], pd.Series]
        The data for the histogram.
    name : Optional[str], optional
        The given name for the plot.
    xbins : Optional[Dict[str, Any]], optional
        Specification of histogram bins (see Plotly documentation for details).
    histnorm : Optional[str], optional
        Specifies the type of normalization for the histogram.
    orientation : str, optional
        The orientation of the histogram ('v' for vertical, 'h' for horizontal).

    Returns
    -------
    go.Histogram
        The Plotly figure object containing the histogram.
    """
    return go.Histogram(
        x=x, xbins=xbins, histnorm=histnorm, orientation=orientation
    )

def generate_plot(
    traces: List[go.Figure], 
    title: Optional[str] = None, 
    xaxis_title: Optional[str] = None, 
    yaxis_title: Optional[str] = None
) -> go.Figure:
    """
    Generate a plot using Plotly based on input traces.

    Parameters
    ----------
    traces : List[go.Figure]
        Contains the list of traces.
    title : Optional[str], optional
        Title for the plot.
    xaxis_title : Optional[str], optional
        x axis title.
    yaxis_title : Optional[str], optional
        y axis title.

    Returns
    -------
    go.Figure
        Returns a generated plot.
    """
    # Create layout
    layout: go.Layout = go.Layout(
        title=title,
        xaxis=dict(title=xaxis_title),
        yaxis=dict(title=yaxis_title)
    )

    # Create figure
    fig: go.Figure = go.Figure(data=traces, layout=layout)
    
    return fig

def write_fig_to_html(fig: go.Figure, output_path: str, filename: str) -> None:
    """
    Writes the fig plotly object to HTML file.

    Parameters
    ----------
    fig : go.Figure
        Figure object to write to HTML file.
    output_path : str
        Path where to save the HTML file.
    filename : str
        Name for the HTML file.

    Returns
    -------
    None
        Writes an HTML file containing the plot to the output directory.
    """
    html_path: str = os.path.join(output_path, f'{filename}.html')
    fig.write_html(html_path)

def plotly_update_layout(
    fig: go.Figure, 
    title: Optional[str] = None, 
    xaxis_title: Optional[str] = None, 
    yaxis_title: Optional[str] = None, 
    legend: Optional[Dict[str, Any]] = None, 
    font: Optional[Dict[str, Any]] = None, 
    margin: Optional[Dict[str, int]] = None, 
    width: Optional[int] = None, 
    height: Optional[int] = None, 
    template: Optional[Union[str, Dict[str, Any]]] = None, 
    plot_bgcolor: Optional[str] = None, 
    paper_bgcolor: Optional[str] = None, 
    annotations: Optional[List[Dict[str, Any]]] = None, 
    shapes: Optional[List[Dict[str, Any]]] = None, 
    images: Optional[List[Dict[str, Any]]] = None, 
    updatemenus: Optional[List[Dict[str, Any]]] = None, 
    sliders: Optional[List[Dict[str, Any]]] = None, 
    scene: Optional[Dict[str, Any]] = None, 
    geo: Optional[Dict[str, Any]] = None, 
    showlegend: Optional[bool] = None
) -> go.Figure:
    """
    Update plotly layout.

    Parameters
    ----------
    fig : go.Figure
        Plotly figure object.
    title : Optional[str], optional
        Title of the plot.
    xaxis_title : Optional[str], optional
        Title of the x-axis.
    yaxis_title : Optional[str], optional
        Title of the y-axis.
    legend : Optional[Dict[str, Any]], optional
        Dictionary specifying legend properties.
    font : Optional[Dict[str, Any]], optional
        Dictionary specifying font properties.
    margin : Optional[Dict[str, int]], optional
        Dictionary specifying margin properties.
    width : Optional[int], optional
        Width of the plot.
    height : Optional[int], optional
        Height of the plot.
    template : Optional[Union[str, Dict[str, Any]]], optional
        Plotly template name or template object to apply to the plot.
    plot_bgcolor : Optional[str], optional
        Background color of the plot.
    paper_bgcolor : Optional[str], optional
        Background color of the plot paper (the area outside the plot).
    annotations : Optional[List[Dict[str, Any]]], optional
        List of dictionaries specifying annotations to be added to the plot.
    shapes : Optional[List[Dict[str, Any]]], optional
        List of dictionaries specifying shapes to be added to the plot.
    images : Optional[List[Dict[str, Any]]], optional
        List of dictionaries specifying images to be added to the plot.
    updatemenus : Optional[List[Dict[str, Any]]], optional
        List of dictionaries specifying update menus to be added to the plot.
    sliders : Optional[List[Dict[str, Any]]], optional
        List of dictionaries specifying sliders to be added to the plot.
    scene : Optional[Dict[str, Any]], optional
        Dictionary specifying the properties of 3D scenes.
    geo : Optional[Dict[str, Any]], optional
        Dictionary specifying the properties of geographic maps.
    showlegend : Optional[bool], optional
        Whether to show the legend.

    Returns
    -------
    go.Figure
        Returns the updated plot.
    """
    return fig.update_layout(
        title=title,
        xaxis_title=xaxis_title,
        yaxis_title=yaxis_title,
        legend=legend,
        font=font,
        margin=margin,
        width=width,
        height=height,
        template=template,
        plot_bgcolor=plot_bgcolor,
        paper_bgcolor=paper_bgcolor,
        annotations=annotations,
        shapes=shapes,
        images=images,
        updatemenus=updatemenus,
        sliders=sliders,
        scene=scene,
        geo=geo,
        showlegend=showlegend
    )

def create_subplots(
    traces: List[go.Figure], 
    rows: int, 
    columns: int, 
    share_x: Optional[bool] = None, 
    share_y: Optional[bool] = None, 
    subplot_titles: Optional[List[str]] = None, 
    coords: Optional[List[List[int]]] = None
) -> go.Figure:
    """
    Creates subplots from traces using Plotly.

    Parameters
    ----------
    traces : List[go.Figure]
        Contains the list of traces.
    rows : int
        Number of rows.
    columns : int
        Number of columns.
    share_x : Optional[bool], optional
        If the traces share the same x-axis.
    share_y : Optional[bool], optional
        If the traces share the same y-axis.
    subplot_titles : Optional[List[str]], optional
        Contains the list of subplot titles.
    coords : Optional[List[List[int]]], optional
        Contains the list of lists of coordinates of rows and columns to add each subplot.

    Returns
    -------
    go.Figure
        Returns a Figure object containing the subplots.
    """
    fig: Figure = make_subplots(
        rows=rows,
        cols=columns,
        shared_xaxes=share_x,
        shared_yaxes=share_y,
        subplot_titles=subplot_titles
    )
    
    for i, trace in enumerate(traces):
        if coords:
            row = coords[i][0]
            column = coords[i][1]
        else:
            row = None
            column = None

        fig.add_trace(trace, row=row, col=column)
    
    return fig

def save_plots_to_html(figures: List[go.Figure], output_path: str, filename: str) -> None:
    """
    Saves a list of plotly.graph_objects objects to one HTML file.

    Parameters
    ----------
    figures : List[go.Figure]
        List that contains the figures to save.
    output_path : str
        Path to the output directory.
    filename : str
        Name for the HTML file.

    Returns
    -------
    None
        Writes an HTML file containing the plots to the output directory.
    """
    html_path: str = os.path.join(output_path, f'{filename}.html')
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