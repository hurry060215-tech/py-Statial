"""Smoke tests — does py-statial import and run at all?"""

import pytest


def test_import():
    import statial
    assert hasattr(statial, "__version__")
    assert statial.__version__ == "0.1.2"


def test_import_functions():
    from statial import (
        get_distances,
        get_abundances,
        calc_contamination,
        Kontextual,
        kontext_curve,
        kontext_plot,
        calc_state_changes,
        make_window,
        parent_combinations,
        get_parent_phylo,
        prep_matrix,
        get_marker_means,
        relabel,
        relabel_kontextual,
        is_kontextual,
    )


def test_make_window_square():
    import numpy as np
    from statial import make_window

    data = {"x": np.array([0, 10, 20]), "y": np.array([0, 10, 20])}
    w = make_window(data, window="square")
    assert w["type"] == "square"
    assert w["xrange"] == (0.0, 20.0)


def test_make_window_convex():
    import numpy as np
    from statial import make_window

    data = {"x": np.array([0, 10, 20, 10]), "y": np.array([0, 20, 0, 10])}
    w = make_window(data, window="convex")
    assert w["type"] == "convex"
    assert "polygon" in w


def test_parent_combinations():
    from statial import parent_combinations

    result = parent_combinations(
        all_types=["A", "B", "C"],
        tcells=["A", "B"],
        tissue=["C"],
    )
    assert len(result) > 0
    assert "from" in result.columns
    assert "to" in result.columns


def test_is_kontextual():
    import pandas as pd
    from statial import is_kontextual

    df = pd.DataFrame({
        "imageID": ["1"],
        "test": ["A__B"],
        "kontextual": [0.5],
        "r": [50],
    })
    assert is_kontextual(df)

    bad = pd.DataFrame({"x": [1]})
    assert not is_kontextual(bad)


def test_relabel():
    import pandas as pd
    import numpy as np
    from statial import relabel

    df = pd.DataFrame({
        "cellType": ["A", "A", "B", "B"],
        "x": [1, 2, 3, 4],
        "y": [1, 2, 3, 4],
    })
    relabeled = relabel(df, labels=["A", "B"], seed=42)
    assert set(relabeled["cellType"]) == {"A", "B"}
    assert len(relabeled) == 4


def test_prep_matrix_kontextual():
    import pandas as pd
    from statial import prep_matrix

    df = pd.DataFrame({
        "imageID": ["1", "1", "2", "2"],
        "test": ["A__B", "A__C", "A__B", "A__C"],
        "kontextual": [0.1, 0.2, 0.3, 0.4],
        "r": [50, 50, 50, 50],
    })
    mat = prep_matrix(df)
    assert mat.shape == (2, 2)


def test_kontextual_basic(small_spatial_data):
    """Basic Kontextual run on synthetic data."""
    from statial import Kontextual

    result = Kontextual(
        cells=small_spatial_data,
        r=100,
        from_types="TypeA",
        to_types="TypeB",
        parent=["TypeA", "TypeB"],
    )
    assert len(result) > 0
    assert "kontextual" in result.columns
    assert "original" in result.columns
