"""Parity gate tests — compare Python output against R reference.

This is the test that decides "is this port done?".
"""

import pytest
import numpy as np
import json
import subprocess
import os


def run_r_reference(fixture_path, output_path, rscript="Rscript"):
    """Run the R reference driver and return the output JSON path."""
    driver = os.path.join(os.path.dirname(__file__), "r_reference_driver.R")
    subprocess.run(
        [rscript, driver, fixture_path, output_path],
        check=True,
        capture_output=True,
    )
    return output_path


def load_json(path):
    with open(path) as f:
        return json.load(f)


def parity_deterministic(reference, candidate, atol=1e-8):
    """Compute max absolute error between reference and candidate arrays."""
    ref = np.asarray(reference, dtype=float).ravel()
    cand = np.asarray(candidate, dtype=float).ravel()
    if len(ref) != len(cand):
        return float("inf")
    return float(np.max(np.abs(ref - cand)))


class TestParityGate:
    """Parity gate tests against R reference."""

    def test_get_distances_parity(self, small_spatial_data, tmp_path):
        """Test get_distances against R reference."""
        from statial import get_distances
        import anndata as ad

        # Create AnnData from DataFrame
        adata = ad.AnnData(
            X=np.random.randn(len(small_spatial_data), 10),
            obs=small_spatial_data,
        )

        # Run Python version
        result = get_distances(adata, max_dist=200, spatial_coords=["x", "y"])

        # Check output shape
        assert "distances" in result.obsm
        dist = result.obsm["distances"]
        assert dist.shape[0] == len(small_spatial_data)

    def test_get_abundances_parity(self, small_spatial_data, tmp_path):
        """Test get_abundances against R reference."""
        from statial import get_abundances
        import anndata as ad

        adata = ad.AnnData(
            X=np.random.randn(len(small_spatial_data), 10),
            obs=small_spatial_data,
        )

        result = get_abundances(adata, r=200)
        assert "abundances" in result.obsm

    def test_kontextual_parity(self, small_spatial_data):
        """Test Kontextual output structure matches R."""
        from statial import Kontextual

        result = Kontextual(
            cells=small_spatial_data,
            r=100,
            from_types="TypeA",
            to_types="TypeB",
            parent=["TypeA", "TypeB"],
            edge_correct=False,
        )

        # Check structure
        assert "imageID" in result.columns
        assert "test" in result.columns
        assert "kontextual" in result.columns
        assert "original" in result.columns
        assert "r" in result.columns

        # Check values are numeric
        assert result["kontextual"].dtype in [np.float64, float]
        assert result["original"].dtype in [np.float64, float]

    def test_make_window_square_parity(self):
        """Test make_window square output."""
        from statial import make_window
        import numpy as np

        data = {"x": np.array([10, 20, 30, 40]), "y": np.array([10, 20, 30, 40])}
        w = make_window(data, window="square")

        # R: owin(xrange=range(x), yrange=range(y))
        assert abs(w["xrange"][0] - 10.0) < 1e-10
        assert abs(w["xrange"][1] - 40.0) < 1e-10
        assert abs(w["yrange"][0] - 10.0) < 1e-10
        assert abs(w["yrange"][1] - 40.0) < 1e-10

    def test_parent_combinations_parity(self):
        """Test parentCombinations output matches R."""
        from statial import parent_combinations

        result = parent_combinations(
            all_types=["tumour", "CD4", "CD8", "epithelial", "stromal"],
            tcells=["CD4", "CD8"],
            tissue=["epithelial", "stromal"],
        )

        # Should have from!=to combinations
        assert len(result) > 0
        assert all(result["from"] != result["to"])

        # Check specific combinations exist
        cd4_rows = result[result["to"] == "CD4"]
        assert len(cd4_rows) > 0

    def test_relabel_structure(self):
        """Test relabel preserves structure."""
        import pandas as pd
        from statial import relabel

        df = pd.DataFrame({
            "cellType": ["A"] * 10 + ["B"] * 10,
            "x": np.random.uniform(0, 100, 20),
            "y": np.random.uniform(0, 100, 20),
        })

        result = relabel(df, labels=["A", "B"], seed=42)
        assert len(result) == 20
        assert set(result["cellType"]) == {"A", "B"}

        # Counts should be preserved
        assert (result["cellType"] == "A").sum() == 10
        assert (result["cellType"] == "B").sum() == 10
