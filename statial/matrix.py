"""Matrix preparation and marker means utilities."""

from __future__ import annotations

import numpy as np
import pandas as pd
import anndata as ad


def prep_matrix(
    result: pd.DataFrame,
    replace_val: float = 0,
    column: str | None = None,
    test: str | None = None,
) -> pd.DataFrame:
    """Convert Kontextual or state changes result to a matrix.

    Mirrors R's ``Statial::prepMatrix``.

    Parameters
    ----------
    result : DataFrame from Kontextual or calcStateChanges
    replace_val : value to replace NAs with
    column : column for values (state changes only)
    test : column for test names (state changes only)

    Returns
    -------
    DataFrame with imageID as index, tests as columns
    """
    if "kontextual" in result.columns:
        mat = result.pivot_table(
            index="imageID",
            columns="test",
            values="kontextual",
            fill_value=replace_val,
        )
        mat.attrs["kontextualResult"] = True
        return mat

    if "primaryCellType" in result.columns:
        if column is not None:
            result = result.copy()
            result["type"] = result[column]
        else:
            result = result.copy()
            result["type"] = result["coef"]

        result["test"] = (
            result["primaryCellType"].astype(str) + "__" +
            result["otherCellType"].astype(str) + "__" +
            result["marker"].astype(str)
        )

        mat = result.pivot_table(
            index="imageID",
            columns="test",
            values="type",
            fill_value=replace_val,
        )
        return mat

    return result


def get_marker_means(
    data: ad.AnnData,
    image_id: str | None = None,
    cell_type: str | None = None,
    region: str | None = None,
    markers: list[str] | None = None,
    assay: int = 0,
    replace_val: float = 0,
) -> pd.DataFrame:
    """Extract average marker expression per cell type per region.

    Mirrors R's ``Statial::getMarkerMeans``.

    Parameters
    ----------
    data : AnnData
    image_id, cell_type, region : column names in .obs
    markers : marker names (default: all var_names)
    assay : which layer
    replace_val : fill value for missing

    Returns
    -------
    DataFrame with means
    """
    if markers is None:
        markers = list(data.var_names)

    cd = data.obs.copy()

    # Get expression
    expr = pd.DataFrame(data.X, index=data.obs_names, columns=data.var_names)[markers]

    use_cols = []
    if image_id is not None:
        cd["imageID"] = cd[image_id]
        use_cols.append("imageID")
    if cell_type is not None:
        cd["cellType"] = cd[cell_type]
        use_cols.append("cellType")
    if region is not None:
        cd["region"] = cd[region]
        use_cols.append("region")

    df = pd.concat([cd[use_cols], expr], axis=1)

    # Melt
    id_vars = use_cols
    value_vars = markers
    melted = df.melt(id_vars=id_vars, value_vars=value_vars, var_name="markers")

    # Group and pivot
    if use_cols:
        grouped = melted.groupby(use_cols + ["markers"])["value"].mean().reset_index()

        if image_id is not None:
            pivot = grouped.pivot_table(
                index="imageID",
                columns=[c for c in use_cols if c != "imageID"] + ["markers"],
                values="value",
                fill_value=replace_val,
            )
            pivot.columns = ["__".join(str(c) for c in col) for col in pivot.columns]
            return pivot
        else:
            pivot = grouped.pivot_table(
                index=use_cols[0] if len(use_cols) == 1 else use_cols,
                columns="markers",
                values="value",
                fill_value=replace_val,
            )
            return pivot

    return expr.mean().to_frame().T
