import matplotlib.pyplot as plt


def _get_edge_verts(reader, tract_id, which_edge="bottom"):
    """Get the vertices of the specified edge of a tract."""
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
    """Plot the outer boundary of a tract."""
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
    """Generate a title for the plot based on the tract and patch IDs."""
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
    for tract_id, patch_id in tract_patch_ids:
        verts = reader.get_patch_vertices(tract_id, patch_id)
        if verts is not None:
            ax.fill(*zip(*verts, strict=True), alpha=0.5, label=f"Tract {tract_id}, Patch {patch_id}")
            ra, dec = zip(*verts, strict=True)
            min_ra, max_ra = min(min_ra, *ra), max(max_ra, *ra)
            min_dec, max_dec = min(min_dec, *dec), max(max_dec, *dec)
    return min_ra, max_ra, min_dec, max_dec


def plot_patches(reader, tract_patch_ids, margin=0.01, tract_outer_boundaries=None, plot_title=None):
    """Plot multiple patches in a single figure."""
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
