"""This module contains example plotting functions for documentation purposes."""

import astropy.units as u
import matplotlib.pyplot as plt
import numpy as np
from astropy.coordinates import SkyCoord

colors = plt.get_cmap("tab10").colors


def demo_rings_plot():
    # Define the number of rings and their latitudes
    n_rings = 5
    lat_edges = np.linspace(-np.pi / 2, np.pi / 2, n_rings + 1)
    lon = np.linspace(-np.pi, np.pi, 300)

    fig = plt.figure(figsize=(8, 5))
    ax = fig.add_subplot(111, projection="mollweide")

    # colors = plt.cm.viridis(np.linspace(0, 1, n_rings))

    for i in range(n_rings):
        lat1 = lat_edges[i]
        lat2 = lat_edges[i + 1]
        # Top edge
        x_top = lon
        y_top = np.full_like(lon, lat2)
        # Bottom edge
        x_bottom = lon[::-1]
        y_bottom = np.full_like(lon, lat1)
        # Combine edges to make a closed polygon
        x_ring = np.concatenate([x_top, x_bottom])
        y_ring = np.concatenate([y_top, y_bottom])
        ax.fill(x_ring, y_ring, color=colors[i], alpha=0.4, edgecolor="black", linewidth=1)

    # Set axis labels and ticks
    ax.set_xlabel("Right Ascension")
    ax.set_ylabel("Declination")
    ax.grid(True, linestyle="--", alpha=0.5)

    plt.tight_layout()
    plt.show()


def demo_rings_tracts_plot():
    # Define the number of rings and their latitudes
    n_rings = 5
    lat_edges = np.linspace(-np.pi / 2, np.pi / 2, n_rings + 1)

    fig = plt.figure(figsize=(8, 5))
    ax = fig.add_subplot(111, projection="mollweide")

    # colors = plt.cm.viridis(np.linspace(0, 1, n_rings))

    for i in range(n_rings):
        lat1 = lat_edges[i]
        lat2 = lat_edges[i + 1]
        # ring_height = lat2 - lat1

        # Define number of quads based on distance from poles
        # Polar rings (0 and 4) get 0 quads, equatorial ring (2) gets most
        if i == 0 or i == 4:  # Polar rings
            n_quads = 0
        elif i == 1 or i == 3:  # Mid-latitude rings
            n_quads = 8
        else:  # Equatorial ring (i == 2)
            n_quads = 14

        if n_quads == 0:
            # Draw the original ring without quads
            lon = np.linspace(-np.pi, np.pi, 300)
            # Top edge
            x_top = lon
            y_top = np.full_like(lon, lat2)
            # Bottom edge
            x_bottom = lon[::-1]
            y_bottom = np.full_like(lon, lat1)
            # Combine edges to make a closed polygon
            x_ring = np.concatenate([x_top, x_bottom])
            y_ring = np.concatenate([y_top, y_bottom])
            ax.fill(x_ring, y_ring, color=colors[i], alpha=0.4, edgecolor="black", linewidth=1)
        else:
            # Create equally spaced longitude edges for the quads
            lon_edges = np.linspace(-np.pi, np.pi, n_quads + 1)

            # Draw each quad in the ring
            for j in range(n_quads):
                lon1 = lon_edges[j]
                lon2 = lon_edges[j + 1]

                # Create quad vertices (rectangular shape)
                x_quad = [lon1, lon2, lon2, lon1, lon1]
                y_quad = [lat1, lat1, lat2, lat2, lat1]

                ax.fill(x_quad, y_quad, color=colors[i], alpha=0.4, edgecolor="black", linewidth=1)

    # Set axis labels and ticks
    ax.set_xlabel("Right Ascension")
    ax.set_ylabel("Declination")
    ax.grid(True, linestyle="--", alpha=0.5)

    plt.tight_layout()
    plt.show()


def _plot_patches(ax, all_patch_verts, single_color=None, alpha=0.5, marker_size=5):
    """Plot the patches on the given axes."""
    legend_added = False  # Flag to ensure we only add one legend entry

    for patch_verts in all_patch_verts:
        ra, dec = zip(*patch_verts, strict=True)
        skycoord = SkyCoord(ra=ra * u.degree, dec=dec * u.degree, frame="icrs")
        ra_deg = skycoord.ra.wrap_at(360 * u.deg).deg
        dec_deg = skycoord.dec.deg

        # Close the polygon
        ra_deg = np.append(ra_deg, ra_deg[0])
        dec_deg = np.append(dec_deg, dec_deg[0])

        # Determine label for legend (only for first patch)
        label = "Patches" if not legend_added else None
        if not legend_added:
            legend_added = True

        # Plot the patches with a single color if specified
        if single_color is not None:
            ax.plot(
                ra_deg,
                dec_deg,
                color=single_color,
                marker="o",
                linestyle="-",
                alpha=alpha,
                markersize=marker_size,
                label=label,
            )
        else:
            ax.plot(
                ra_deg,
                dec_deg,
                marker="o",
                linestyle="-",
                alpha=alpha,
                markersize=marker_size,
                label=label,
            )


def _plot_tract(ax, tract_id, tract_verts, color):
    """Plot the tract boundary."""
    ra, dec = zip(*tract_verts, strict=True)
    skycoord = SkyCoord(ra=ra * u.degree, dec=dec * u.degree, frame="icrs")
    ra_deg = skycoord.ra.wrap_at(360 * u.deg).deg
    dec_deg = skycoord.dec.deg

    # Close the polygon
    ra_deg = np.append(ra_deg, ra_deg[0])
    dec_deg = np.append(dec_deg, dec_deg[0])

    # Plot the tract boundary
    ax.fill(
        ra_deg,
        dec_deg,
        color=color,
        alpha=0.3,
        edgecolor="black",
        linewidth=1.5,
        label=f"Tract {tract_id} (Inner)",
    )


def _get_ra_dec_range(patches):
    """Get the RA/Dec range from the patches."""
    min_ra, max_ra = float("inf"), float("-inf")
    min_dec, max_dec = float("inf"), float("-inf")
    for patch in patches:
        ra, dec = zip(*patch, strict=True)
        min_ra = min(min(ra), min_ra)
        max_ra = max(max(ra), max_ra)
        min_dec = min(min(dec), min_dec)
        max_dec = max(max(dec), max_dec)
    return (min_ra, max_ra, min_dec, max_dec)


def _plot_data(ax, data=None, target_tract=None):
    """Plot the data points on the given axes."""
    if data is not None and target_tract is not None:
        for i in range(len(data)):
            row = data.iloc[i]
            ra, dec = row["ra"], row["dec"]
            color = colors[9] if row["tract"] == target_tract else colors[1]
            ax.plot(ra, dec, marker=".", color=color, markersize=5, alpha=0.7)


def plot_patches_in_tract(reader, tract_id, data=None):
    """Plot the patches in a specific tract with zoomed-in view."""
    fig, ax = plt.subplots(figsize=(10, 6))

    # Plot the patches and tract boundaries
    all_patch_verts = []
    for patch_id in range(100):
        patch_verts = reader.get_patch_vertices(tract_id, patch_id)
        if patch_verts is None:
            raise ValueError(f"Patch {patch_id} not found in tract {tract_id}.")
        all_patch_verts.append(patch_verts)
    _plot_patches(ax, all_patch_verts, single_color=colors[0])
    tract_verts = reader.get_tract_vertices(tract_id)
    _plot_tract(ax, tract_id, tract_verts, color=colors[7])

    # Plot the data, if provided
    _plot_data(ax, data, target_tract=tract_id)

    # Set zoom level based on the RA/Dec range of the patches
    min_ra, max_ra, min_dec, max_dec = _get_ra_dec_range(all_patch_verts)
    ax.set_xlim(min_ra - 1, max_ra + 1)
    ax.set_ylim(min_dec - 0.25, max_dec + 0.25)

    # Add a legend for the tract boundaries, labels, and title
    ax.legend(loc="upper right", fontsize="small")
    ax.set_xlabel("RA (deg)")
    ax.set_ylabel("Dec (deg)")
    ax.set_title(f"Zoomed-In View of DP1 Subset Tract Designation for Tract {tract_id}")
    ax.grid(True)

    ax.legend()
    plt.tight_layout()
    plt.show()
