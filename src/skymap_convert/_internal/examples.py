"""This module contains example plotting functions for documentation purposes."""

import matplotlib.pyplot as plt
import numpy as np


def demo_rings_plot():
    # Define the number of rings and their latitudes
    n_rings = 5
    lat_edges = np.linspace(-np.pi / 2, np.pi / 2, n_rings + 1)
    lon = np.linspace(-np.pi, np.pi, 300)

    fig = plt.figure(figsize=(8, 5))
    ax = fig.add_subplot(111, projection="mollweide")

    colors = plt.cm.viridis(np.linspace(0, 1, n_rings))

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

    colors = plt.cm.viridis(np.linspace(0, 1, n_rings))

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
