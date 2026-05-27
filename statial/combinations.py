"""Create parent-child cell type combinations."""

from __future__ import annotations

import pandas as pd
from itertools import product


def parent_combinations(
    all_types: list[str],
    parent_list: dict[str, list[str]] | None = None,
    **kwargs,
) -> pd.DataFrame:
    """Create all pairwise cell type relationships from parent definitions.

    Mirrors R's ``Statial::parentCombinations``.

    Parameters
    ----------
    all_types : list of all cell types
    parent_list : dict mapping parent name → list of children
    **kwargs : alternative way to pass parent vectors (name=vector)

    Returns
    -------
    DataFrame with columns: from, to, parent, parent_name
    """
    if parent_list is None:
        parent_list = kwargs

    rows = []
    for parent_name, children in parent_list.items():
        for child in children:
            for from_type in all_types:
                if from_type != child:
                    rows.append({
                        "from": from_type,
                        "to": child,
                        "parent": children,
                        "parent_name": parent_name,
                    })

    return pd.DataFrame(rows, columns=["from", "to", "parent", "parent_name"])


def get_parent_phylo(phylo_tree) -> dict[str, list[str]]:
    """Extract parent-children relationships from a phylogenetic tree.

    Mirrors R's ``Statial::getParentPhylo``.

    Parameters
    ----------
    phylo_tree : dict with 'tip.label' and 'edge' keys (phylo format)

    Returns
    -------
    dict mapping parent name → list of children
    """
    if isinstance(phylo_tree, dict) and "clust_tree" in phylo_tree:
        phylo_tree = phylo_tree["clust_tree"]

    tip_labels = phylo_tree["tip.label"]
    edge = phylo_tree["edge"]

    n_tips = len(tip_labels)
    node_labels = {i + 1: tip_labels[i] for i in range(n_tips)}

    # Find internal nodes with >1 child (among tips)
    parent_children = {}
    for parent_node, child_node in edge:
        if child_node <= n_tips:
            child_name = node_labels[child_node]
            if parent_node not in parent_children:
                parent_children[parent_node] = []
            parent_children[parent_node].append(child_name)

    # Filter to parents with >1 child
    result = {}
    for i, (node, children) in enumerate(parent_children.items()):
        if len(children) > 1:
            result[f"parent_{i+1}"] = children

    return result
