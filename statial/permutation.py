"""Cell label permutation for significance testing."""

from __future__ import annotations

import numpy as np
import pandas as pd


def relabel(
    image: pd.DataFrame,
    labels: list[str] | None = None,
    seed: int | None = None,
) -> pd.DataFrame:
    """Permute cell type labels within an image.

    Mirrors R's ``Statial::relabel``.

    Parameters
    ----------
    image : DataFrame with cellType column
    labels : cell types to permute (None = all)
    seed : random seed

    Returns
    -------
    DataFrame with permuted cellType
    """
    if seed is not None:
        np.random.seed(seed)

    result = image.copy()

    if labels is None:
        result["cellType"] = np.random.permutation(result["cellType"].values)
    else:
        mask = result["cellType"].isin(labels)
        to_relabel = result[mask].copy()
        not_relabel = result[~mask]

        to_relabel["cellType"] = np.random.permutation(to_relabel["cellType"].values)
        result = pd.concat([to_relabel, not_relabel])

    return result


def relabel_kontextual(
    cells,
    n_sim: int = 1,
    r=None,
    from_types=None,
    to_types=None,
    parent=None,
    image=None,
    return_images: bool = False,
    inhom: bool = True,
    edge: bool = False,
    cores: int = 1,
    spatial_coords=None,
    cell_type: str = "cellType",
    image_id: str = "imageID",
    **kwargs,
) -> pd.DataFrame:
    """Permute labels and calculate Kontextual values.

    Mirrors R's ``Statial::relabelKontextual``.

    Parameters
    ----------
    cells : DataFrame or AnnData
    n_sim : number of permutations
    r : radius or list of radii
    from_types, to_types, parent : cell type specifications
    image : subset to specific images
    return_images : return permuted images too
    inhom, edge : Kontextual parameters

    Returns
    -------
    DataFrame with Kontextual values for original and permuted data
    """
    from .kontextual import Kontextual

    if spatial_coords is None:
        spatial_coords = ["x", "y"]

    # Normalize input
    if isinstance(cells, pd.DataFrame):
        df = cells.copy()
        if "cellID" not in df.columns:
            df["cellID"] = range(len(df))
    else:
        raise ValueError("cells must be a DataFrame")

    if r is None:
        r = np.arange(10, 110, 10)

    parent_df = pd.DataFrame({
        "from": [from_types],
        "to": [to_types],
        "parent": [parent],
    })

    # Original
    orig = Kontextual(
        cells=df,
        r=r,
        parent_df=parent_df,
        image=image,
        inhom=inhom,
        edge_correct=edge,
        include_original=True,
        spatial_coords=spatial_coords,
        cell_type=cell_type,
        image_id=image_id,
    )
    orig["type"] = "original"
    orig["imageID"] = 1

    results = [orig]

    # Permutations
    for sim in range(n_sim):
        perm_df = relabel(df, labels=parent, seed=42 + sim)
        perm_result = Kontextual(
            cells=perm_df,
            r=r,
            parent_df=parent_df,
            image=image,
            inhom=inhom,
            edge_correct=edge,
            include_original=True,
            spatial_coords=spatial_coords,
            cell_type=cell_type,
            image_id=image_id,
        )
        perm_result["type"] = "randomised"
        perm_result["imageID"] = str(sim + 2)
        results.append(perm_result)

    all_results = pd.concat(results, ignore_index=True)

    if return_images:
        return all_results, None

    return all_results
