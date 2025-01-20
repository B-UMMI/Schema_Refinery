import os
import pandas as pd
import plotly.graph_objects as go
import plotly.offline
from plotly.subplots import make_subplots
from typing import List, Optional, Union, Dict, Any

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
    fig: go.Figure = make_subplots(
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
