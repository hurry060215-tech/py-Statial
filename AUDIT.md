# AUDIT.md - R Function Coverage

## Exported functions

| R function | Python function | Status |
|---|---|---|
| `Kontextual` | `Kontextual` | ported |
| `calcContamination` | `calc_contamination` | ported |
| `calcStateChanges` | `calc_state_changes` | ported |
| `getAbundances` | `get_abundances` | ported |
| `getDistances` | `get_distances` | ported |
| `getMarkerMeans` | `get_marker_means` | ported |
| `getParentPhylo` | `get_parent_phylo` | ported |
| `isKontextual` | `is_kontextual` | ported |
| `kontextCurve` | `kontext_curve` | ported |
| `kontextPlot` | `kontext_plot` | ported |
| `makeWindow` | `make_window` | ported |
| `parentCombinations` | `parent_combinations` | ported |
| `plotStateChanges` | `---` | MISSING |
| `prepMatrix` | `prep_matrix` | ported |
| `relabel` | `relabel` | ported |
| `relabelKontextual` | `relabel_kontextual` | ported |

**Total: 16 exported, 15 ported, 1 missing**

## Internal (non-exported) functions

| R function | Python function | Status |
|---|---|---|
| `preProcessing` | `---` | missing-or-inlined |
| `validateDf` | `---` | missing-or-inlined |

## Notes

- `plotStateChanges`: Visualization-only function (ggplot2), excluded from parity gate
- `kontextPlot`: Ported as matplotlib-based function
- All core computational functions are ported with parity validation