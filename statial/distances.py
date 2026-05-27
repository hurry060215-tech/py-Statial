"""Calculate pairwise distances between cell types."""

from __future__ import annotations

import numpy as np
import pandas as pd
import anndata as ad

from .utils import nearest_distances


def get_distances(
    cells: ad.AnnData,
    max_dist: float | None = None,
    image_id: str = "imageID",
    spatial_coords: list[str] | None = None,
    cell_type: str = "cellType",
    red_dim_name: str = "distances",
    dist_fun: str = "min",
) -> ad.AnnData:
    """Calculate euclidean distance from each cell to nearest cell of each type.

    Mirrors R's ``Statial::getDistances``.

    Parameters
    ----------
    cells : AnnData with spatial coordinates in .obs
    max_dist : maximum distance considered (default: max of coordinate range)
    image_id : column in .obs with image identifiers
    spatial_coords : columns with x, y coordinates (default: ["x", "y"])
    cell_type : column in .obs with cell type labels
    red_dim_name : key in .obsm to store results
    dist_fun : "min" or "abundance"

    Returns
    -------
    AnnData with distances stored in .obsm[red_dim_name]
    """
    if spatial_coords is None:
        spatial_coords = ["x", "y"]

    cd = cells.obs.copy()
    ct_col = cd[cell_type].values
    img_col = cd[image_id].values if image_id in cd.columns else np.zeros(len(cd), dtype=int)
    coords = cd[spatial_coords].values.astype(float)
    cell_ids = cd.index.values

    if max_dist is None:
        max_dist = max(
            coords[:, 0].max() - coords[:, 0].min(),
            coords[:, 1].max() - coords[:, 1].min(),
        )

    unique_types = sorted(cd[cell_type].unique())
    unique_images = sorted(pd.unique(img_col))

    all_dfs = []
    for img in unique_images:
        mask = img_col == img
        img_coords = coords[mask]
        img_labels = ct_col[mask]
        img_ids = cell_ids[mask]

        df = nearest_distances(img_coords, img_labels, unique_types, max_dist, dist_fun)
        df.index = img_ids
        all_dfs.append(df)

    distances = pd.concat(all_dfs)
    distances = distances.loc[cell_ids]

    cells.obsm[red_dim_name] = distances.values
    cells.uns[f"{red_dim_name}_columns"] = list(distances.columns)

    return cells
