# Discovery — py-statial

## Direct check

**No existing omicverse port found.** Safe to start a new port.

Checked: `github.com/omicverse` for `py-statial`, `py-Statial`, `statial` — 0 results.

## Dependency audit

### R dependencies (from DESCRIPTION)

| R dep | omicverse mirror | Decision |
|---|---|---|
| `BiocParallel` | none | Native Python (multiprocessing) |
| `spatstat.geom` | none | Native Python (scipy.spatial.cKDTree) |
| `concaveman` | none | Native Python (scipy.spatial.ConvexHull fallback) |
| `data.table` | none | Native Python (pandas) |
| `spatstat.explore` | none | Not needed (used for exploration only) |
| `dplyr` | none | Native Python (pandas) |
| `tidyr` | none | Native Python (pandas) |
| `SingleCellExperiment` | none | Native Python (anndata.AnnData) |
| `tibble` | none | Native Python (pandas) |
| `stringr` | none | Native Python (str methods) |
| `tidyselect` | none | Native Python (pandas column selection) |
| `ggplot2` | none | matplotlib (for kontextPlot) |
| `plotly` | none | Skipped (interactive plotting) |
| `purrr` | none | Native Python (itertools/functools) |
| `ranger` | none | Native Python (sklearn RandomForestClassifier) |
| `limma` | none | Native Python (numpy.linalg.lstsq) |
| `SpatialExperiment` | none | Native Python (anndata.AnnData) |
| `cluster` | none | Not needed (R-specific clustering) |
| `treekoR` | none | Not needed (used only by getParentPhylo) |
| `edgeR` | none | Native Python (scipy.stats) |
| `S4Vectors` | none | Native Python (pandas) |
| `SummarizedExperiment` | none | Native Python (anndata.AnnData) |

### Suggests (optional)

| R dep | Decision |
|---|---|
| `BiocStyle` | Out of scope (R formatting) |
| `knitr` | Out of scope (R vignettes) |
| `testthat` | Replaced by pytest |
| `ClassifyR` | Not needed |
| `spicyR` | Not needed |
| `ggsurvfit` | Not needed |
| `lisaClust` | Not needed |
| `survival` | Not needed |

## Decisions

- All R dependencies replaced by native Python equivalents (scipy, sklearn, pandas, anndata)
- No omicverse mirrors found or needed
- `plotStateChanges` skipped (visualization-only function)
- `concaveman` replaced by convex hull (concave window uses convex as fallback)
