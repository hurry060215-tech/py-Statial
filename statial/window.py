"""Create observation windows for spatial point patterns."""

from __future__ import annotations

import numpy as np
from scipy.spatial import ConvexHull


def make_window(
    data: dict | None = None,
    x: np.ndarray | None = None,
    y: np.ndarray | None = None,
    window: str = "square",
    window_length: float | None = None,
) -> dict:
    """Create a window definition for spatial data.

    Mirrors R's ``Statial::makeWindow``.

    Parameters
    ----------
    data : dict with 'x' and 'y' keys, or DataFrame-like
    x, y : coordinate arrays (used if data is None)
    window : "square", "convex", or "concave"
    window_length : tuning parameter for concave windows

    Returns
    -------
    dict with window type and polygon/bounds info
    """
    if data is not None:
        x = np.asarray(data["x"]) if isinstance(data, dict) else np.asarray(data.iloc[:, 0])
        y = np.asarray(data["y"]) if isinstance(data, dict) else np.asarray(data.iloc[:, 1])
    if x is None or y is None:
        raise ValueError("Must provide data or x,y coordinates")

    x_range = (float(np.min(x)), float(np.max(x)))
    y_range = (float(np.min(y)), float(np.max(y)))

    if window == "square":
        return {"type": "square", "xrange": x_range, "yrange": y_range}

    if window == "convex":
        points = np.column_stack([x, y])
        hull = ConvexHull(points)
        polygon = points[hull.vertices]
        return {"type": "convex", "polygon": polygon, "xrange": x_range, "yrange": y_range}

    if window == "concave":
        # Simplified concave window using alpha shape concept
        if window_length is None or np.isnan(window_length):
            wl = (x_range[1] - x_range[0]) / 20.0
        else:
            wl = (x_range[1] - x_range[0]) / 20.0 * window_length

        points = np.column_stack([x, y])
        # Use convex hull as fallback (concaveman is R-specific)
        hull = ConvexHull(points)
        polygon = points[hull.vertices]
        return {
            "type": "concave",
            "polygon": polygon,
            "window_length": wl,
            "xrange": x_range,
            "yrange": y_range,
        }

    raise ValueError(f"Unknown window type: {window}")
