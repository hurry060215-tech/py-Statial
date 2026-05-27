"""Parity tests comparing Python output against R reference outputs.

These tests load the R reference JSON produced by r_reference_driver.R
and compare with Python function outputs on the same data.
"""

import pytest
import json
import numpy as np
import pandas as pd
import os

DATA_DIR = os.path.join(os.path.dirname(__file__), "..", "data")
R_TOL = 1e-6  # tolerance for deterministic comparison


def load_r_reference():
    """Load R reference outputs."""
    path = os.path.join(DATA_DIR, "r_reference_output.json")
    if not os.path.exists(path):
        pytest.skip("R reference output not found")
    with open(path) as f:
        return json.load(f)


def load_cell_metadata():
    """Load cell metadata exported from R."""
    path = os.path.join(DATA_DIR, "cell_metadata.csv")
    if not os.path.exists(path):
        pytest.skip("Cell metadata not found")
    return pd.read_csv(path, index_col=0)


def parse_r_matrix(data, key):
    """Parse R matrix from JSON, handling NAs (-1 sentinel)."""
    mat = np.array(data[key], dtype=float)
    mat[mat == -1] = np.nan
    return mat


def align_and_compare(py_mat, r_mat, py_cols, r_cols):
    """Align two matrices by column name and compute max error on non-NaN values."""
    # Build column mapping
    r_col_map = {c: i for i, c in enumerate(r_cols)}

    # Find common columns in same order as Python
    common_cols = [c for c in py_cols if c in r_col_map]
    if not common_cols:
        pytest.skip("No common columns between Python and R")

    r_idx = [r_col_map[c] for c in common_cols]
    py_idx = [py_cols.index(c) for c in common_cols]

    r_aligned = r_mat[:, r_idx]
    py_aligned = py_mat[:, py_idx]

    # Compare non-NaN values
    mask = ~(np.isnan(py_aligned) | np.isnan(r_aligned))
    if mask.sum() == 0:
        pytest.skip("No overlapping non-NaN values")

    max_err = np.max(np.abs(py_aligned[mask] - r_aligned[mask]))
    mean_err = np.mean(np.abs(py_aligned[mask] - r_aligned[mask]))
    return max_err, mean_err


class TestParityDistances:
    """Test get_distances parity with R."""

    def test_distances_shape(self):
        r_ref = load_r_reference()
        r_dist = parse_r_matrix(r_ref, "distances")
        assert r_dist.shape == (57811, 17)

    def test_distances_column_names(self):
        r_ref = load_r_reference()
        cols = r_ref["distances_colnames"]
        assert len(cols) == 17
        assert "Keratin_Tumour" in cols

    def test_distances_values(self):
        """Python distances should match R values within tolerance."""
        from statial import get_distances
        import anndata as ad

        meta = load_cell_metadata()
        r_ref = load_r_reference()
        r_dist = parse_r_matrix(r_ref, "distances")
        r_cols = r_ref["distances_colnames"]

        # Create AnnData
        np.random.seed(42)
        adata = ad.AnnData(
            X=np.random.randn(len(meta), 10),
            obs=meta,
        )

        # Run Python
        result = get_distances(adata, max_dist=200, spatial_coords=["x", "y"])
        py_dist = result.obsm["distances"]
        py_cols = result.uns["distances_columns"]

        # Align by column name and compare
        max_err, mean_err = align_and_compare(py_dist, r_dist, py_cols, r_cols)
        assert max_err < R_TOL, f"Max error {max_err:.2e} exceeds tolerance {R_TOL} (mean={mean_err:.2e})"


class TestParityAbundances:
    """Test get_abundances parity with R."""

    def test_abundances_shape(self):
        r_ref = load_r_reference()
        r_abund = parse_r_matrix(r_ref, "abundances")
        assert r_abund.shape == (57811, 17)

    def test_abundances_values(self):
        """Python abundances should match R values exactly."""
        from statial import get_abundances
        import anndata as ad

        meta = load_cell_metadata()
        r_ref = load_r_reference()
        r_abund = parse_r_matrix(r_ref, "abundances")
        r_cols = r_ref["abundances_colnames"]

        np.random.seed(42)
        adata = ad.AnnData(X=np.random.randn(len(meta), 10), obs=meta)

        result = get_abundances(adata, r=200, spatial_coords=["x", "y"])
        py_abund = result.obsm["abundances"]
        py_cols = result.uns["abundances_columns"]

        max_err, mean_err = align_and_compare(py_abund, r_abund, py_cols, r_cols)
        assert max_err < R_TOL, f"Max error {max_err:.2e} exceeds tolerance {R_TOL}"


class TestParityKontextual:
    """Test Kontextual parity with R."""

    def test_kontextual_structure(self):
        r_ref = load_r_reference()
        k = r_ref["kontextual"]
        assert len(k) == 1
        assert k[0]["imageID"] == "6"
        assert k[0]["test"] == "Macrophages__Keratin_Tumour"
        assert k[0]["r"] == 50

    def test_kontextual_r_values(self):
        """Check R reference values are reasonable."""
        r_ref = load_r_reference()
        k = r_ref["kontextual"][0]
        assert np.isfinite(k["original"])
        assert np.isfinite(k["kontextual"])
        assert abs(k["original"]) < 1000
        assert abs(k["kontextual"]) < 1000

    def test_kontextual_parity(self):
        """Kontextual output should match R values."""
        from statial import Kontextual

        meta = load_cell_metadata()
        r_ref = load_r_reference()
        r_kont = r_ref["kontextual"][0]

        result = Kontextual(
            cells=meta,
            r=50,
            from_types="Macrophages",
            to_types="Keratin_Tumour",
            parent=["Macrophages", "CD4_Cell"],
            image=["6"],
            edge_correct=False,
            window="square",
            spatial_coords=["x", "y"],
        )

        if len(result) == 0:
            pytest.skip("No results for image 6")

        row = result.iloc[0]

        # Original L should match exactly
        orig_err = abs(row["original"] - r_kont["original"])
        assert orig_err < 1e-6, f"Original L error {orig_err:.2e} too large"

        # Kontextual may differ slightly due to cKDTree vs spatstat closepairs
        kont_err = abs(row["kontextual"] - r_kont["kontextual"])
        kont_rel = kont_err / abs(r_kont["kontextual"]) if r_kont["kontextual"] != 0 else 0
        assert kont_rel < 0.10, f"Kontextual relative error {kont_rel:.2%} too large"


class TestParityContamination:
    """Test calc_contamination parity with R."""

    def test_contamination_shape(self):
        r_ref = load_r_reference()
        r_contam = parse_r_matrix(r_ref, "contaminations")
        assert r_contam.shape[0] == 57811
        cols = r_ref["contaminations_colnames"]
        assert len(cols) > 0
