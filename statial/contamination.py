"""Calculate marker contamination scores using random forest."""

from __future__ import annotations

import numpy as np
import pandas as pd
import anndata as ad
from sklearn.ensemble import RandomForestClassifier


def calc_contamination(
    cells: ad.AnnData,
    markers: list[str] | None = None,
    num_trees: int = 100,
    verbose: bool = False,
    missing_replacement: float = 0,
    assay: str = "intensities",
    cell_type: str = "cellType",
    red_dim_name: str = "contaminations",
) -> ad.AnnData:
    """Calculate contamination scores using random forest classification.

    Mirrors R's ``Statial::calcContamination``.

    Parameters
    ----------
    cells : AnnData with marker expression in .X or .layers[assay]
    markers : list of marker names (default: all)
    num_trees : number of trees in random forest
    verbose : print model info
    missing_replacement : value for missing markers
    assay : which layer to use for expression
    cell_type : column in .obs with cell types
    red_dim_name : key in .obsm for output

    Returns
    -------
    AnnData with contamination scores in .obsm[red_dim_name]
    """
    if markers is None:
        if hasattr(cells, "var_names"):
            markers = list(cells.var_names)
        else:
            markers = list(range(cells.shape[1]))

    # Get expression data
    if assay in cells.layers:
        expr = pd.DataFrame(cells.layers[assay], index=cells.obs_names, columns=cells.var_names)
    else:
        expr = pd.DataFrame(cells.X, index=cells.obs_names, columns=cells.var_names)

    # Subset to markers that exist
    available_markers = [m for m in markers if m in expr.columns]
    if not available_markers:
        available_markers = list(expr.columns)

    expr = expr[available_markers].copy()

    # Replace NaN/missing
    expr = expr.fillna(missing_replacement)
    expr = expr.replace([np.inf, -np.inf], missing_replacement)

    # Prepare training data
    y = cells.obs[cell_type].values
    rf_data = expr.copy()

    # Fit random forest
    clf = RandomForestClassifier(
        n_estimators=num_trees,
        oob_score=True,
        random_state=42,
        n_jobs=-1,
    )
    clf.fit(rf_data, y)

    if verbose:
        print(f"OOB score: {clf.oob_score_:.4f}")

    # Get class probabilities
    predictions = clf.predict_proba(rf_data)
    class_names = clf.classes_

    # Build output DataFrame
    pred_df = pd.DataFrame(predictions, index=cells.obs_names, columns=class_names)

    # Max probability
    pred_df["rfMaxCellProb"] = pred_df.max(axis=1)

    # Second largest probability
    sorted_probs = np.sort(predictions, axis=1)
    pred_df["rfSecondLargestCellProb"] = sorted_probs[:, -2]

    # Main cell probability (probability assigned to true class)
    main_probs = np.zeros(len(y))
    for i, ct in enumerate(y):
        if ct in class_names:
            idx = list(class_names).index(ct)
            main_probs[i] = predictions[i, idx]
    pred_df["rfMainCellProb"] = main_probs

    # Store result
    cells.obsm[red_dim_name] = pred_df

    return cells
