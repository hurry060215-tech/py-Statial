"""Validation utilities."""

from __future__ import annotations

import pandas as pd


def is_kontextual(obj) -> bool:
    """Test whether an object is a kontextualResult.

    Mirrors R's ``Statial::isKontextual``.

    Parameters
    ----------
    obj : object to test

    Returns
    -------
    True if obj has the expected kontextualResult columns
    """
    if not isinstance(obj, pd.DataFrame):
        return False

    required = {"imageID", "test", "kontextual", "r"}
    return required.issubset(set(obj.columns))
