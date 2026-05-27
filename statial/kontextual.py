"""Kontextual: conditional spatial relationship analysis."""

from __future__ import annotations

import numpy as np
import pandas as pd
from itertools import product


def _l_function(close_pairs_df, counts_df, child1, child2, r, area):
    """Calculate centered L function.

    Mirrors R's ``Statial:::.Lfunction``.
    """
    n_child1 = (counts_df["cellTypeI"] == child1).sum()
    n_child2 = (counts_df["cellTypeI"] == child2).sum()

    numerator = counts_df.loc[counts_df["cellTypeI"] == child1, child2].sum()

    lambda2 = n_child2 / area
    if lambda2 == 0 or n_child1 == 0:
        return np.nan

    kvalue = (numerator / lambda2) * (1 / n_child1)
    centered_l = np.sqrt(kvalue / np.pi) - r

    return centered_l


def _l_inhom_function(close_pairs_df, counts_df, child1, child2, r, area):
    """Calculate inhomogeneous centered L function.

    Mirrors R's ``Statial:::.Linhomfunction``.
    """
    n_child1 = (counts_df["cellTypeI"] == child1).sum()

    weight_child1 = counts_df[child1].values / (np.pi * r**2)
    weight_child2 = counts_df[child2].values / (np.pi * r**2)

    mask = (close_pairs_df["cellTypeI"] == child1) & (close_pairs_df["cellTypeJ"] == child2)
    cp = close_pairs_df[mask].copy()

    if len(cp) == 0:
        return np.nan

    # Map i to index in counts_df for weight lookup
    idx_map = {v: i for i, v in enumerate(counts_df["i"].values)}
    wi = np.array([weight_child1[idx_map.get(ii, 0)] for ii in cp["i"].values])
    wj = np.array([weight_child2[idx_map.get(jj, 0)] for jj in cp["j"].values])

    cp_vals = cp["edge"].values / (wi * wj)
    numerator = cp_vals.sum()

    child1_mask = counts_df["cellTypeI"] == child1
    denominator = (1.0 / weight_child1[child1_mask]).sum()

    if denominator == 0:
        return np.nan

    centered_l = np.sqrt((numerator / denominator) / np.pi) - r
    return centered_l


def _kontext(close_pairs_df, counts_df, child1, child2, parent, r, return_weight=False):
    """Calculate Kontextual value.

    Mirrors R's ``Statial:::.Kontext``.

    Parameters
    ----------
    close_pairs_df : DataFrame of close pairs (already filtered by radius)
    counts_df : DataFrame of counts per cell (already filtered by radius). If None, built from close_pairs_df.
    """
    if counts_df is None:
        # Build counts from the (already filtered) close pairs
        agg = close_pairs_df.groupby(["i", "cellTypeI", "cellTypeJ"]).agg(n=("edge", "sum")).reset_index()
        counts_df = agg.pivot_table(index=["i", "cellTypeI"], columns="cellTypeJ", values="n", fill_value=0).reset_index()
        counts_df.columns.name = None

    n_parent = counts_df["cellTypeI"].isin(parent).sum()
    n_child1 = (counts_df["cellTypeI"] == child1).sum()
    n_child2 = (counts_df["cellTypeI"] == child2).sum()

    # Parent columns that exist in counts
    parent_names = [p for p in parent if p in counts_df.columns]
    if not parent_names:
        return np.nan

    counts_parent = counts_df[parent_names].sum(axis=1).values
    lambda_parent = counts_parent / (np.pi * r**2)
    lambda_child1 = lambda_parent.copy()
    lambda_child2 = (lambda_parent / n_parent) * n_child2 if n_parent > 0 else np.zeros_like(lambda_parent)

    mask = (close_pairs_df["cellTypeI"] == child1) & (close_pairs_df["cellTypeJ"] == child2)
    cp = close_pairs_df[mask].copy()

    if len(cp) == 0:
        return np.nan

    idx_map = {v: i for i, v in enumerate(counts_df["i"].values)}
    li = np.array([lambda_child1[idx_map.get(ii, 0)] for ii in cp["i"].values])
    lj = np.array([lambda_child2[idx_map.get(jj, 0)] for jj in cp["j"].values])

    cp["weightParent"] = (cp["edge"].values * li) / lj
    numerator = cp["weightParent"].sum()

    child1_mask = counts_df["cellTypeI"] == child1
    denominator = lambda_child1[child1_mask].sum()

    if denominator == 0:
        return np.nan

    if return_weight:
        return {"lambdaChild1": lambda_child1, "lambdaChild2": lambda_child2, "denominator": denominator}

    centered_l = np.sqrt(numerator / denominator / np.pi) - r
    return centered_l


def _border_edge(coords, window_info, max_d):
    """Edge correction for border cells.

    Simplified version of R's ``Statial:::.borderEdge``.
    """
    n = len(coords)
    edge = np.ones(n)

    x_range = window_info["xrange"]
    y_range = window_info["yrange"]

    for i in range(n):
        x, y = coords[i]
        # Check if cell is within max_d of any border
        dist_to_border = min(
            x - x_range[0],
            x_range[1] - x,
            y - y_range[0],
            y_range[1] - y,
        )
        if dist_to_border < max_d:
            # Approximate edge correction: fraction of circle inside window
            # Simplified: linear interpolation
            frac = max(0.0, min(1.0, dist_to_border / max_d))
            edge[i] = frac if frac > 0 else 0.01  # avoid division by zero

    return edge


def Kontextual(
    cells,
    r,
    parent_df=None,
    from_types=None,
    to_types=None,
    parent=None,
    image=None,
    inhom=False,
    edge_correct=True,
    window="convex",
    window_length=None,
    include_original=True,
    spatial_coords=None,
    cell_type="cellType",
    image_id="imageID",
    cores=1,
):
    """Evaluate pairwise cell relationships, conditional on a 3rd population.

    Mirrors R's ``Statial::Kontextual``.

    Parameters
    ----------
    cells : DataFrame or dict of DataFrames with cellType, imageID, x, y columns
    r : radius or list of radii
    parent_df : DataFrame with from, to, parent columns
    from_types, to_types, parent : alternative to parent_df
    image : subset to specific images
    inhom : use inhomogeneous L function
    edge_correct : perform edge correction
    window : "square", "convex", or "concave"
    include_original : include original L values in output
    spatial_coords : column names for x, y
    cell_type, image_id : column names
    cores : number of cores (unused in Python, kept for API compat)

    Returns
    -------
    DataFrame with imageID, test, original, kontextual, r columns
    """
    if spatial_coords is None:
        spatial_coords = ["x", "y"]

    # Build parent_df from from/to/parent args if needed
    if parent_df is None:
        if from_types is None or to_types is None or parent is None:
            raise ValueError("Must specify parent_df or from_types/to_types/parent")
        parent_df = pd.DataFrame({
            "from": [from_types],
            "to": [to_types],
            "parent": [parent],
        })

    # Normalize input to list of DataFrames
    if isinstance(cells, pd.DataFrame):
        df = cells.copy()
        df["cellID"] = range(len(df))
        df = df.rename(columns={
            spatial_coords[0]: "x",
            spatial_coords[1]: "y",
            cell_type: "cellType",
            image_id: "imageID",
        })
        # Ensure consistent types for imageID comparison
        df["imageID"] = df["imageID"].astype(str)
        if image is not None:
            image_str = [str(i) for i in image]
            df = df[df["imageID"].isin(image_str)]
        images = {img: sub for img, sub in df.groupby("imageID")}
    elif isinstance(cells, dict):
        images = cells
    else:
        raise ValueError("cells must be DataFrame or dict of DataFrames")

    r_vals = np.sort(np.atleast_1d(r))

    # Build image info
    results = []
    for img_name, img_df in images.items():
        img_df = img_df.copy()
        img_df["cellID"] = range(len(img_df))
        coords = img_df[["x", "y"]].values
        labels = img_df["cellType"].values
        cell_ids = img_df["cellID"].values

        # Compute window
        win = make_window_simple(coords, window, window_length)
        area = (win["xrange"][1] - win["xrange"][0]) * (win["yrange"][1] - win["yrange"][0])

        # Find close pairs for max radius
        max_r = r_vals[-1]
        cp = _find_close_pairs(coords, labels, cell_ids, max_r)

        if len(cp) == 0:
            for _, row in parent_df.iterrows():
                for rv in r_vals:
                    results.append({
                        "imageID": img_name,
                        "test": f"{row['from']}__{row['to']}",
                        "original": np.nan,
                        "kontextual": np.nan,
                        "r": rv,
                        "inhomL": inhom,
                    })
            continue

        # Edge correction
        if edge_correct:
            edge = _border_edge(coords, win, max_r)
            edge_map = {cell_ids[i]: edge[i] for i in range(len(cell_ids))}
            cp["edge"] = cp["i"].map(edge_map)
        else:
            cp["edge"] = 1.0

        # Process each radius
        for rv in r_vals:
            cp_r = cp[cp["d"] < rv].copy() if rv < max_r else cp.copy()

            # Build counts
            if len(cp_r) > 0:
                agg = cp_r.groupby(["i", "cellTypeI", "cellTypeJ"]).agg(n=("edge", "sum")).reset_index()
                counts = agg.pivot_table(index=["i", "cellTypeI"], columns="cellTypeJ", values="n", fill_value=0).reset_index()
                counts.columns.name = None
            else:
                counts = pd.DataFrame(columns=["i", "cellTypeI"])

            # Process each parent combination
            for _, prow in parent_df.iterrows():
                from_t = prow["from"]
                to_t = prow["to"]
                par = prow["parent"] if isinstance(prow["parent"], list) else [prow["parent"]]

                # Check if cell types exist
                unique_types = set(labels)
                if from_t not in unique_types or to_t not in unique_types:
                    results.append({
                        "imageID": img_name,
                        "test": f"{from_t}__{to_t}",
                        "original": np.nan,
                        "kontextual": np.nan,
                        "r": rv,
                        "inhomL": inhom,
                    })
                    continue

                # Original L function
                if inhom:
                    orig = _l_inhom_function(cp_r, counts, from_t, to_t, rv, area)
                else:
                    orig = _l_function(cp_r, counts, from_t, to_t, rv, area)

                # Kontextual value
                kont = _kontext(cp_r, counts, from_t, to_t, par, rv)

                results.append({
                    "imageID": img_name,
                    "test": f"{from_t}__{to_t}",
                    "original": orig,
                    "kontextual": kont,
                    "r": rv,
                    "inhomL": inhom,
                })

    result_df = pd.DataFrame(results)

    if not include_original and len(result_df) > 0:
        result_df = result_df[["imageID", "test", "kontextual", "r"]]

    return result_df


def _find_close_pairs(coords, labels, cell_ids, r_max):
    """Find all cell pairs within r_max distance."""
    from scipy.spatial import cKDTree

    if len(coords) < 2:
        return pd.DataFrame(columns=["i", "j", "d", "cellTypeI", "cellTypeJ", "edge"])

    tree = cKDTree(coords)
    pairs = tree.query_pairs(r_max, output_type="ndarray")

    if len(pairs) == 0:
        return pd.DataFrame(columns=["i", "j", "d", "cellTypeI", "cellTypeJ", "edge"])

    i_idx = pairs[:, 0]
    j_idx = pairs[:, 1]
    d = np.linalg.norm(coords[i_idx] - coords[j_idx], axis=1)

    # Both directions
    all_i = np.concatenate([i_idx, j_idx])
    all_j = np.concatenate([j_idx, i_idx])
    all_d = np.concatenate([d, d])

    return pd.DataFrame({
        "i": cell_ids[all_i],
        "j": cell_ids[all_j],
        "d": all_d,
        "cellTypeI": labels[all_i],
        "cellTypeJ": labels[all_j],
        "edge": 1.0,  # placeholder, set later
    })


def make_window_simple(coords, window="square", window_length=None):
    """Simplified window creation."""
    x_range = (float(coords[:, 0].min()), float(coords[:, 0].max()))
    y_range = (float(coords[:, 1].min()), float(coords[:, 1].max()))
    return {"type": window, "xrange": x_range, "yrange": y_range}


def kontext_curve(
    cells,
    from_types,
    to_types,
    parent,
    image=None,
    rs=None,
    inhom=False,
    edge=True,
    se=False,
    n_sim=20,
    cores=1,
    image_id="imageID",
    cell_type="cellType",
    **kwargs,
):
    """Evaluate Kontextual over a range of radii.

    Mirrors R's ``Statial::kontextCurve``.
    """
    if rs is None:
        rs = np.arange(10, 110, 10)

    result = Kontextual(
        cells=cells,
        r=rs,
        from_types=from_types,
        to_types=to_types,
        parent=parent,
        image=image,
        inhom=inhom,
        edge_correct=edge,
        cores=cores,
        image_id=image_id,
        cell_type=cell_type,
        include_original=True,
        **kwargs,
    )

    if result.empty:
        return result

    rs_df = result[["r", "original", "kontextual"]].copy()

    if se:
        from .permutation import relabel_kontextual
        se_df = relabel_kontextual(
            cells=cells,
            n_sim=n_sim,
            r=rs,
            from_types=from_types,
            to_types=to_types,
            parent=parent,
            image=image,
            inhom=inhom,
            edge=edge,
            cores=cores,
            image_id=image_id,
            cell_type=cell_type,
        )
        se_df = se_df[se_df["type"] != "original"]
        se_stats = se_df.groupby("r").agg(
            originalSd=("original", "sd"),
            kontextualSd=("kontextual", "sd"),
        ).reset_index()
        rs_df = rs_df.merge(se_stats, on="r", how="left")

    return rs_df


def kontext_plot(rs_df):
    """Plot Kontextual values over radii.

    Mirrors R's ``Statial::kontextPlot``.

    Returns matplotlib figure.
    """
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(figsize=(8, 5))

    has_sd = "kontextualSd" in rs_df.columns

    if has_sd:
        ax.fill_between(
            rs_df["r"],
            rs_df["kontextual"] - rs_df["kontextualSd"],
            rs_df["kontextual"] + rs_df["kontextualSd"],
            alpha=0.2, color="blue",
        )
        ax.fill_between(
            rs_df["r"],
            rs_df["original"] - rs_df["originalSd"],
            rs_df["original"] + rs_df["originalSd"],
            alpha=0.2, color="orange",
        )

    ax.plot(rs_df["r"], rs_df["kontextual"], "o-", label="Kontextual", color="blue")
    ax.plot(rs_df["r"], rs_df["original"], "o-", label="Original L", color="orange")
    ax.axhline(y=0, linestyle="--", color="red", alpha=0.5)
    ax.set_xlabel("Radius (r)")
    ax.set_ylabel("Relationship value")
    ax.legend()
    ax.set_title("Kontextual vs Original L values")

    return fig
