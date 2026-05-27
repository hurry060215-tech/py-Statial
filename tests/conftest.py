"""Pytest fixtures for py-statial tests."""

import pytest
import numpy as np
import pandas as pd


@pytest.fixture
def small_spatial_data():
    """Create a small synthetic spatial dataset for testing."""
    np.random.seed(42)
    n = 200

    x = np.random.uniform(0, 1000, n)
    y = np.random.uniform(0, 1000, n)

    cell_types = np.random.choice(
        ["TypeA", "TypeB", "TypeC"],
        size=n,
        p=[0.4, 0.35, 0.25],
    )

    df = pd.DataFrame({
        "x": x,
        "y": y,
        "cellType": cell_types,
        "imageID": "image1",
    })
    df.index = [f"cell_{i}" for i in range(n)]

    return df


@pytest.fixture
def multi_image_data():
    """Create multi-image spatial dataset."""
    np.random.seed(42)
    dfs = []
    for img in ["img1", "img2"]:
        n = 100
        x = np.random.uniform(0, 500, n)
        y = np.random.uniform(0, 500, n)
        ct = np.random.choice(["A", "B", "C"], size=n)
        df = pd.DataFrame({"x": x, "y": y, "cellType": ct, "imageID": img})
        df.index = [f"{img}_cell_{i}" for i in range(n)]
        dfs.append(df)
    return pd.concat(dfs)


@pytest.fixture
def keren_sce_path():
    """Path to kerenSCE fixture if available."""
    import os
    path = os.path.join(os.path.dirname(__file__), "..", "data", "kerenSCE.rds")
    if os.path.exists(path):
        return path
    return None
