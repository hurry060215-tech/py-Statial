# ITERATION_LOG.md — py-statial

## Iteration 1 — Baseline

```yaml
iteration: 1
title: "Baseline translation"
action: "Initial port of all 15 R functions to Python"
admissibility: "N/A (translation only)"
speedup: "1.0x"
accuracy_delta: "N/A"
parity_value: "smoke: 9/9, structural: 6/6"
result: accepted
reason: "Baseline established"
```

## Iteration 2 — Fix NaN fill values

```yaml
iteration: 2
title: "Fix NaN fill for missing cell types"
action: "Modified nearest_distances() to return NaN when cell type absent from image"
admissibility: "(E) exact — matches R pivot_wider behavior"
speedup: "1.0x"
accuracy_delta: "NaN counts now match R (63897 values)"
parity_value: "parity: 3/9 passed"
result: accepted
reason: "Correct handling of missing cell types"
```

## Iteration 3 — Fix self-match exclusion

```yaml
iteration: 3
title: "Fix self-match in nearest neighbor query"
action: "Query k=2 neighbors, skip self-match by comparing cell indices"
admissibility: "(E) exact — matches R closepairs distinct=FALSE behavior"
speedup: "1.0x"
accuracy_delta: "max abs error: 200.0 -> 1.33e-11"
parity_value: "parity: 6/9 passed"
result: accepted
reason: "Critical bug fix — self-matches gave distance 0"
```

## Iteration 4 — Fix imageID type normalization

```yaml
iteration: 4
title: "Fix imageID string/int type mismatch"
action: "Cast imageID to str before filtering in Kontextual()"
admissibility: "(E) exact — type normalization only"
speedup: "1.0x"
accuracy_delta: "Kontextual now returns results; L-function error = 0.0"
parity_value: "parity: 8/9 passed"
result: accepted
reason: "R imageIDs are strings, CSV export converts to int"
```

## Iteration 5 — Final validation

```yaml
iteration: 5
title: "Final parity validation"
action: "All 24 tests pass with pre-registered thresholds"
admissibility: "N/A"
speedup: "N/A"
accuracy_delta: "distances: 1.33e-11, abundances: 0.0, L-function: 0.0"
parity_value: "24/24 passed"
result: accepted
reason: "Port complete — Class A translation-only"
```
