"""Shared utility functions for py-statial."""

from __future__ import annotations

import numpy as np
import pandas as pd
from scipy.spatial import cKDTree


def nearest_distances(
    coords: np.ndarray,
    labels: np.ndarray,
    cell_types: list[str],
    max_dist: float,
    dist_fun: str = "min",
) -> pd.DataFrame:
    """For each cell, compute distance to nearest cell of each type.

    Parameters
    ----------
    coords : (n, 2) array of x, y coordinates
    labels : (n,) array of cell type labels
    cell_types : list of unique cell types
    max_dist : maximum distance to consider
    dist_fun : "min" for nearest distance, "abundance" for count within max_dist

    Returns
    -------
    DataFrame with one column per cell type, rows aligned to input order.
    NaN for cell types not present in the image.
    """
    n = len(coords)
    tree = cKDTree(coords)

    result = {}
    for ct in cell_types:
        mask = labels == ct
        if not mask.any():
            result[ct] = np.full(n, np.nan)
            continue

        ct_indices = np.where(mask)[0]
        ct_tree = cKDTree(coords[mask])

        if dist_fun == "min":
            # Query k=2 to handle self-matches vectorized
            dists, idxs = ct_tree.query(coords, k=2, distance_upper_bound=max_dist)
            # Clamp invalid indices (inf returns n_trees as index)
            idxs_safe = np.clip(idxs, 0, len(ct_indices) - 1)
            real_idxs = ct_indices[idxs_safe]
            self_match = real_idxs == np.arange(n)[:, None]
            # Invalid entries (inf distance) are not self-matches
            self_match[dists == np.inf] = False
            # Take first non-self neighbor; if first is self, take second
            nearest = np.where(self_match[:, 0], dists[:, 1], dists[:, 0])
            # Fill inf with max_dist
            nearest[nearest == np.inf] = max_dist
            result[ct] = nearest

        elif dist_fun == "abundance":
            # Bulk computation via sparse distance matrix
            sp = tree.sparse_distance_matrix(ct_tree, max_dist, output_type="coo_matrix")
            counts = np.array(sp.getnnz(axis=1)).flatten().astype(float)
            # Subtract self-matches
            counts[mask] -= 1
            counts = np.maximum(counts, 0)
            result[ct] = counts

    return pd.DataFrame(result)
