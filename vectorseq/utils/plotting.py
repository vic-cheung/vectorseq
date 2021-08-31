import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns


def surface_plot_3d(
    df: pd.DataFrame,
    title: str = "",
    xlabel: str = "",
    ylabel: str = "",
    zlabel: str = "",
    figsize: tuple = (10, 10),
    elevation: int = None,
    azimuth: int = None,
    cmap: str = "magma",
    edgecolor: str = "white",
    linewidths: float = 0.5,
):
    """
    Create 3D surface plot.  Datapoints are at the vertices (not faces) of the surface.

    Args:
        df: dataframe of dimension (y,x).  columns will become X-axis, index will become Y-axis.
        title: plot title
        (x,y,z)label: labels for x,y,z axis
        figsize: figure size (width, height)
        elevation: up-down rotation of 3D plot
        azimuth: left-right rotation of 3D plot
        cmap: matplotlib color map theme
        edgecolor: color of edges.
        linewidths: width of lines.
    Returns:
        matplotlib figure object
    """
    fig = plt.figure(figsize=figsize)
    ax = plt.axes(projection="3d")
    X = np.tile(df.columns.to_numpy(), (df.shape[0], 1)).astype(float)
    Y = np.tile(df.index.to_numpy(), (df.shape[1], 1)).T.astype(float)
    Z = df.to_numpy()
    ax.plot_surface(
        X,
        Y,
        Z,
        rstride=1,
        cstride=1,
        cmap=cmap,
        edgecolor=edgecolor,
        linewidths=linewidths,
    )
    ax.set_title(title)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_zlabel(zlabel)
    ax.view_init(elevation, azimuth)
    return fig


def heatmap_plot(
    df: pd.DataFrame,
    title: str = "",
    xlabel: str = "",
    ylabel: str = "",
    figsize: tuple = (10, 10),
    cmap: str = "magma",
    edgecolor: str = "white",
    linewidths: float = 0.5,
    annotate: bool = True,
    annotate_format: str = "d",
):
    """
    Create 2D heatmap.  Datapoints are the square faces of the grid.

    Args:
        df: dataframe of dimension (y,x).  columns will become X-axis, index will become Y-axis.
        title: plot title
        (x,y)label: labels for x,y,z axis
        figsize: figure size (width, height)
        cmap: matplotlib color map theme
        edgecolor: color of edges.
        linewidths: width of lines.
        annotate: whether to annotate heatmap faces with numeric value
        annotate_format: string format settings for numeric annotation
            (e.g. "d" for integers, ".2f" for 2 decimals,
            ".4E" for 4-decimal exponential notation)

    Returns:
        matplotlib figure object
    """
    fig = plt.figure(figsize=figsize)
    ax = sns.heatmap(
        data=df,
        annot=annotate,
        fmt=annotate_format,
        cmap=cmap,
        edgecolor=edgecolor,
        linewidths=linewidths,
    )
    ax.set_title(title)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    return fig
