"""Core class for py-statial — wraps AnnData and provides method chaining."""

from __future__ import annotations

import anndata as ad


class Statial:
    """Main class owning an AnnData object with method chaining."""

    def __init__(self, adata: ad.AnnData):
        self.adata = adata

    def get_distances(self, **kwargs) -> "Statial":
        from .distances import get_distances
        self.adata = get_distances(self.adata, **kwargs)
        return self

    def get_abundances(self, **kwargs) -> "Statial":
        from .abundances import get_abundances
        self.adata = get_abundances(self.adata, **kwargs)
        return self

    def calc_contamination(self, **kwargs) -> "Statial":
        from .contamination import calc_contamination
        self.adata = calc_contamination(self.adata, **kwargs)
        return self
