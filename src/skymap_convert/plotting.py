import matplotlib.pyplot as plt


def _get_edge_verts(reader, tract_id, which_edge="bottom"):
    """Get the vertices of the specified edge of a tract.

    Parameters
    ----------
    reader : ConvertedSkymapReader
        Reader instance for accessing skymap data.
    tract_id : int
        ID of the tract to extract edge vertices from.
    which_edge : str, optional
        Which edge of the tract to extract. Must be one of:
        'bottom', 'left', 'top', 'right' (default: 'bottom').

    Returns
    -------
    list of list of float
        List of [RA, Dec] vertex coordinates in degrees for the specified edge.

    Raises
    ------
    ValueError
        If which_edge is not one of the valid options.

    Notes
    -----
    Edge vertices are extracted from patch boundaries assuming a 10x10 patch grid.
    Corner patches are included in both adjacent edges.
    """
    # The bottom edge is defined by the patches: 0, 1, ..., 8, 9
    # We can take just the bottom-right vertex of each, which will be the 1st vertex of each patch.
    # The next edge is the left edge, and if we grab the bottom-left vertex of each patch,
    # this will cover that bottom-left patch's bottom-left vertex.
    # The trick here is each patch_ids list is double-counting the corner patches, ie, the bottom-left
    # patch is included in both the bottom and left edges.
    if which_edge == "bottom":
        patch_ids = range(10)  # Patches 0 to 9
        vertex_index = 1  # Bottom-right vertex of each patch
    elif which_edge == "left":
        patch_ids = range(9, 100, 10)  # Patches 9, 19, 29, ..., 99
        vertex_index = 0  # Bottom-left vertex of each patch
    elif which_edge == "top":
        patch_ids = range(99, 89, -1)  # Patches 99 to 90
        vertex_index = 3  # Top-left vertex of each patch
    elif which_edge == "right":
        patch_ids = range(90, 0, -10)  # Patches 90, 80, ..., 10, 0
        vertex_index = 2  # Top-right vertex of each patch
    else:
        raise ValueError("Invalid edge specified. Choose from 'bottom', 'left', 'top', or 'right'.")

    # Get the vertices for the specified edge
    edge_verts = []
    for i in patch_ids:
        edge_verts.append(reader.get_patch_vertices(tract_id, i)[vertex_index])
    return edge_verts


def _plot_tract_outer_boundary(reader, ax, tract_id, min_ra, max_ra, min_dec, max_dec):
    """Plot the outer boundary of a tract and update coordinate bounds.

    Parameters
    ----------
    reader : ConvertedSkymapReader
        Reader instance for accessing skymap data.
    ax : matplotlib.axes.Axes
        Matplotlib axes object to plot on.
    tract_id : int
        ID of the tract to plot.
    min_ra : float
        Current minimum RA value in degrees.
    max_ra : float
        Current maximum RA value in degrees.
    min_dec : float
        Current minimum Dec value in degrees.
    max_dec : float
        Current maximum Dec value in degrees.

    Returns
    -------
    tuple of float
        Updated (min_ra, max_ra, min_dec, max_dec) bounds in degrees

    Notes
    -----
    Plots the tract boundary as a light gray filled polygon with dashed outline.
    """
    # Get the vertices of the outer boundary of the tract
    tract_outer_verts = []
    bottom_edge = _get_edge_verts(reader, tract_id, which_edge="bottom")
    tract_outer_verts += bottom_edge  # Bottom edge vertices
    left_edge = _get_edge_verts(reader, tract_id, which_edge="left")
    tract_outer_verts += left_edge  # Left edge vertices
    top_edge = _get_edge_verts(reader, tract_id, which_edge="top")
    tract_outer_verts += top_edge  # Top edge vertices
    right_edge = _get_edge_verts(reader, tract_id, which_edge="right")
    tract_outer_verts += right_edge  # Right edge vertices

    # Update the overall min/max values
    ra, dec = zip(*tract_outer_verts, strict=True)
    min_ra, max_ra = min(min_ra, *ra), max(max_ra, *ra)
    min_dec, max_dec = min(min_dec, *dec), max(max_dec, *dec)

    # Draw the outer boundary of the tract
    ax.fill(
        *zip(*tract_outer_verts, strict=True),
        facecolor="lightgray",
        alpha=0.5,
        linestyle="--",
        linewidth=1,
        edgecolor="black",
        label="Tract Outer Region",
    )

    return min_ra, max_ra, min_dec, max_dec


def _generate_title(tract_patch_ids, tract_outer_boundaries):
    """Generate a descriptive title for the plot based on tract and patch IDs.

    Parameters
    ----------
    tract_patch_ids : list of tuple
        List of (tract_id, patch_id) tuples for patches being plotted.
    tract_outer_boundaries : list of int or None
        List of tract IDs whose outer boundaries are being plotted.

    Returns
    -------
    str
        Generated plot title describing the content being displayed.

    Notes
    -----
    Title format varies based on the number of tracts:
    - Single tract: "Skymap Patch Plot for Tract {id}"
    - Multiple tracts: "Skymap Patch Plot for Tracts: {id1}, {id2}, ..."
    """
    if not tract_patch_ids:
        return "No patches to plot"

    tract_ids = set(tract_id for tract_id, _ in tract_patch_ids)
    tract_bound_ids = set(tract_outer_boundaries) if tract_outer_boundaries else set()
    all_tract_ids = tract_ids.union(tract_bound_ids)

    title = "Skymap Patch Plot"
    if len(all_tract_ids) == 1:
        title += f" for Tract {next(iter(all_tract_ids))}"
    elif all_tract_ids:
        title += f" for Tracts: {', '.join(map(str, sorted(all_tract_ids)))}"

    return title


def _plot_patches(reader, ax, tract_patch_ids, min_ra, max_ra, min_dec, max_dec):
    """Plot individual patches and update coordinate bounds.

    Parameters
    ----------
    reader : ConvertedSkymapReader
        Reader instance for accessing skymap data.
    ax : matplotlib.axes.Axes
        Matplotlib axes object to plot on.
    tract_patch_ids : list of tuple
        List of (tract_id, patch_id) tuples specifying patches to plot.
    min_ra : float
        Current minimum RA value in degrees.
    max_ra : float
        Current maximum RA value in degrees.
    min_dec : float
        Current minimum Dec value in degrees.
    max_dec : float
        Current maximum Dec value in degrees.

    Returns
    -------
    tuple of float
        Updated (min_ra, max_ra, min_dec, max_dec) bounds in degrees.

    Notes
    -----
    Each patch is plotted as a filled polygon with 50% transparency.
    """
    for tract_id, patch_id in tract_patch_ids:
        verts = reader.get_patch_vertices(tract_id, patch_id)
        if verts is not None:
            ax.fill(*zip(*verts, strict=True), alpha=0.5, label=f"Tract {tract_id}, Patch {patch_id}")
            ra, dec = zip(*verts, strict=True)
            min_ra, max_ra = min(min_ra, *ra), max(max_ra, *ra)
            min_dec, max_dec = min(min_dec, *dec), max(max_dec, *dec)
    return min_ra, max_ra, min_dec, max_dec


def plot_patches(reader, tract_patch_ids, margin=0.01, tract_outer_boundaries=None, plot_title=None):
    """Plot multiple patches in a single figure with optional tract boundaries.

    Creates a matplotlib figure showing specified skymap patches as filled polygons,
    with optional tract outer boundaries for context.

    Parameters
    ----------
    reader : ConvertedSkymapReader
        Reader instance for accessing skymap data.
    tract_patch_ids : list of tuple
        List of (tract_id, patch_id) tuples specifying patches to plot.
    margin : float, optional
        Margin as a fraction of the plot range to add around the plotted area
        (default: 0.01).
    tract_outer_boundaries : int or list of int, optional
        Tract ID(s) whose outer boundaries should be plotted for context.
        If int, plots boundary for single tract. If list, plots boundaries
        for all specified tracts (default: None).
    plot_title : str, optional
        Custom title for the plot. If None, generates automatic title based
        on plotted content (default: None).

    Returns
    -------
    None
        Displays the plot using matplotlib.pyplot.show().

    Examples
    --------
    Plot patches from a single tract:
    >>> plot_patches(reader, [(0, 45), (0, 46), (0, 55)])

    Plot patches with tract boundary context:
    >>> plot_patches(reader, [(0, 45), (1, 20)], tract_outer_boundaries=[0, 1])

    Plot with custom title and margin:
    >>> plot_patches(reader, [(0, 45)], margin=0.05, plot_title="Custom Title")

    Notes
    -----
    - Patches are plotted as filled polygons with 50% transparency
    - Tract boundaries are shown as light gray dashed outlines
    - Coordinate axes are labeled in degrees (RA/Dec)
    - Legend shows tract and patch IDs for each plotted element
    """
    fig, ax = plt.subplots(figsize=(10, 6))

    # Initialize min/max values for RA and Dec
    min_ra, max_ra = float("inf"), float("-inf")
    min_dec, max_dec = float("inf"), float("-inf")

    # Optionally plot the outer boundaries of the tracts
    if tract_outer_boundaries:
        if not isinstance(tract_outer_boundaries, list):
            tract_outer_boundaries = [tract_outer_boundaries]
        for tract_id in tract_outer_boundaries:
            min_ra, max_ra, min_dec, max_dec = _plot_tract_outer_boundary(
                reader, ax, tract_id, min_ra, max_ra, min_dec, max_dec
            )

    # Iterate over tract-patch pairs and plot each patch
    min_ra, max_ra, min_dec, max_dec = _plot_patches(
        reader, ax, tract_patch_ids, min_ra, max_ra, min_dec, max_dec
    )

    # Set bounds with optional margins
    if margin > 0:
        # Calculate margins based on the min/max RA and Dec
        ra_interval = max_ra - min_ra
        dec_interval = max_dec - min_dec

        margin_ra_amount = ra_interval * margin
        margin_dec_amount = dec_interval * margin

        min_ra -= margin_ra_amount
        max_ra += margin_ra_amount
        min_dec -= margin_dec_amount
        max_dec += margin_dec_amount

    # Set the axis limits and labels
    ax.set_xlim(min_ra, max_ra)
    ax.set_ylim(min_dec, max_dec)
    ax.set_xlabel("RA (degrees)")
    ax.set_ylabel("Dec (degrees)")
    plot_title = plot_title if plot_title else _generate_title(tract_patch_ids, tract_outer_boundaries)
    ax.set_title(plot_title)
    ax.legend()

    # Show the plot
    plt.show()
