"""Calculate cell state changes using linear models."""

from __future__ import annotations

import numpy as np
import pandas as pd
import anndata as ad
from scipy import stats


def calc_state_changes(
    cells: ad.AnnData,
    markers: list[str] | None = None,
    from_types: list[str] | None = None,
    to_types: list[str] | None = None,
    image: list[str] | None = None,
    type_key: str = "distances",
    assay: int = 0,
    cell_type: str = "cellType",
    image_id: str = "imageID",
    contamination: str | None = None,
    test: str = "g",
    min_cells: int = 20,
    verbose: bool = False,
    timeout: int = 10,
) -> pd.DataFrame:
    """Build linear models measuring marker-based state changes.

    Mirrors R's ``Statial::calcStateChanges``.

    Parameters
    ----------
    cells : AnnData with distances in .obsm[type_key]
    markers : list of markers (default: all)
    from_types : primary cell types (default: all)
    to_types : interacting cell types (default: all)
    image : subset to specific images
    type_key : key in .obsm for distances
    assay : which assay/layer index
    cell_type, image_id : column names
    contamination : key for contamination scores in .obsm
    test : "g" for Gaussian, "nb" for negative binomial
    min_cells : minimum cells per group

    Returns
    -------
    DataFrame with columns: imageID, primaryCellType, otherCellType, marker, coef, tval, pval, fdr
    """
    if markers is None:
        markers = list(cells.var_names)

    cd = cells.obs.copy()
    img_col = cd[image_id].values if image_id in cd.columns else np.zeros(len(cd), dtype=int)

    if image is not None:
        mask = np.isin(img_col, image)
        cells = cells[mask]
        cd = cells.obs.copy()
        img_col = cd[image_id].values

    ct_col = cd[cell_type].values

    if to_types is None:
        to_types = sorted(np.unique(ct_col))
    if from_types is None:
        from_types = sorted(np.unique(ct_col))

    # Filter to from_types
    from_mask = np.isin(ct_col, from_types)
    cells = cells[from_mask]
    cd = cells.obs.copy()
    ct_col = cd[cell_type].values
    img_col = cd[image_id].values

    # Get distances
    if type_key in cells.obsm:
        distances = pd.DataFrame(
            cells.obsm[type_key],
            index=cells.obs_names,
        )
        # Subset to to_types columns if available
        if hasattr(cells.uns, "get") and f"{type_key}_columns" in cells.uns:
            dist_cols = cells.uns[f"{type_key}_columns"]
            distances.columns = dist_cols
        if to_types and all(t in distances.columns for t in to_types):
            distances = distances[to_types]
    else:
        raise ValueError(f"Distances not found in .obsm['{type_key}']")

    # Get expression
    if hasattr(cells, "layers") and len(cells.layers) > assay:
        intensities = pd.DataFrame(
            cells.layers[assay] if isinstance(assay, str) else list(cells.layers.values())[assay],
            index=cells.obs_names,
            columns=cells.var_names,
        )
    else:
        intensities = pd.DataFrame(cells.X, index=cells.obs_names, columns=cells.var_names)

    intensities = intensities[markers]

    # Filter by min cells per group
    group_df = pd.DataFrame({"imageID": img_col, "cellType": ct_col})
    group_counts = group_df.groupby(["imageID", "cellType"]).size().reset_index(name="n")
    valid_groups = group_counts[group_counts["n"] > min_cells][["imageID", "cellType"]]

    valid_mask = np.zeros(len(cd), dtype=bool)
    for _, row in valid_groups.iterrows():
        valid_mask |= (img_col == row["imageID"]) & (ct_col == row["cellType"])

    if not valid_mask.any():
        return pd.DataFrame(columns=["imageID", "primaryCellType", "otherCellType", "marker", "coef", "tval", "pval", "fdr"])

    cells = cells[valid_mask]
    cd = cells.obs.copy()
    ct_col = cd[cell_type].values
    img_col = cd[image_id].values

    # Rebuild data
    if type_key in cells.obsm:
        distances = pd.DataFrame(cells.obsm[type_key], index=cells.obs_names)
        if hasattr(cells.uns, "get") and f"{type_key}_columns" in cells.uns:
            distances.columns = cells.uns[f"{type_key}_columns"]
        if to_types and all(t in distances.columns for t in to_types):
            distances = distances[to_types]

    intensities = pd.DataFrame(cells.X, index=cells.obs_names, columns=cells.var_names)[markers]

    # Get contamination scores if specified
    contam_df = None
    if contamination is not None and contamination in cells.obsm:
        contam_df = pd.DataFrame(cells.obsm[contamination], index=cells.obs_names)

    # Split by image+cellType and fit models
    all_results = []
    groups = pd.DataFrame({"imageID": img_col, "cellType": ct_col})

    for (img, ct), idx in groups.groupby(["imageID", "cellType"]).groups.items():
        if len(idx) < min_cells:
            continue

        d_sub = distances.loc[idx]
        i_sub = intensities.loc[idx]

        for other_ct in to_types:
            if other_ct not in d_sub.columns:
                continue

            x = d_sub[other_ct].values
            if len(np.unique(x)) <= 1:
                continue

            for marker in markers:
                y = i_sub[marker].values

                # Build design matrix
                design = np.column_stack([np.ones(len(x)), x])

                # Add contamination covariates if available
                if contam_df is not None:
                    c_sub = contam_df.loc[idx]
                    # Drop all-zero or constant columns
                    c_sub = c_sub.loc[:, c_sub.std() > 0]
                    if len(c_sub.columns) > 0:
                        design = np.column_stack([design, c_sub.values])

                # OLS fit
                try:
                    beta, residuals, rank, sv = np.linalg.lstsq(design, y, rcond=None)
                    y_hat = design @ beta
                    mse = np.sum((y - y_hat) ** 2) / (len(y) - rank)
                    if mse > 0:
                        se = np.sqrt(np.diag(np.linalg.pinv(design.T @ design)) * mse)
                        tval = beta[1] / se[1] if se[1] > 0 else 0
                        pval = 2 * stats.t.sf(abs(tval), df=len(y) - rank)
                    else:
                        tval = 0
                        pval = 1.0

                    all_results.append({
                        "imageID": img,
                        "primaryCellType": ct,
                        "otherCellType": other_ct,
                        "marker": marker,
                        "coef": beta[1],
                        "tval": tval,
                        "pval": pval,
                    })
                except Exception:
                    continue

    if not all_results:
        return pd.DataFrame(columns=["imageID", "primaryCellType", "otherCellType", "marker", "coef", "tval", "pval", "fdr"])

    result = pd.DataFrame(all_results)
    result["fdr"] = stats.false_discovery_control(result["pval"].values, method="bh")
    result = result.sort_values("pval").reset_index(drop=True)

    return result
