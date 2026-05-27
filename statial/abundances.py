"""Calculate cell type abundances (K-function based)."""

from __future__ import annotations

import numpy as np
import pandas as pd
import anndata as ad

from .utils import nearest_distances


def get_abundances(
    cells: ad.AnnData,
    r: float = 200,
    dist_fun: str = "abundance",
    red_dim_name: str = "abundances",
    cell_type: str = "cellType",
    image_id: str = "imageID",
    spatial_coords: list[str] | None = None,
) -> ad.AnnData:
    """Calculate abundance of each cell type around each cell.

    Mirrors R's ``Statial::getAbundances``.

    Parameters
    ----------
    cells : AnnData
    r : radius for abundance calculation
    dist_fun : distance function (default "abundance")
    red_dim_name : key in .obsm
    cell_type : column in .obs
    image_id : column in .obs
    spatial_coords : columns with x, y

    Returns
    -------
    AnnData with abundances in .obsm[red_dim_name]
    """
    if spatial_coords is None:
        spatial_coords = ["x", "y"]

    cd = cells.obs.copy()
    ct_col = cd[cell_type].values
    img_col = cd[image_id].values if image_id in cd.columns else np.zeros(len(cd), dtype=int)
    coords = cd[spatial_coords].values.astype(float)
    cell_ids = cd.index.values

    unique_types = sorted(cd[cell_type].unique())
    unique_images = sorted(pd.unique(img_col))

    all_dfs = []
    for img in unique_images:
        mask = img_col == img
        img_coords = coords[mask]
        img_labels = ct_col[mask]
        img_ids = cell_ids[mask]

        df = nearest_distances(img_coords, img_labels, unique_types, r, dist_fun)
        df.index = img_ids
        all_dfs.append(df)

    abundances = pd.concat(all_dfs)
    abundances = abundances.loc[cell_ids]

    cells.obsm[red_dim_name] = abundances.values
    cells.uns[f"{red_dim_name}_columns"] = list(abundances.columns)

    return cells
