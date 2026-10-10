"""Discover target-preserving aggregate partitions from sparse interactions."""
from .rational import q
from .finite import partition


def aggregate_partition(operator, initial_labels):
    """Coarsest refinement with exact closed group sums.

    Columns with equal aggregate contributions to every current group stay
    together. This is aggregate closure, not a claim of graph isomorphism.
    """
    if len(initial_labels) != operator.n:
        raise ValueError("One initial label per coordinate required")
    classes = partition(initial_labels)
    rounds = []
    while True:
        groups = max(classes)+1
        contributions = [[q(0)]*groups for _ in range(operator.n)]
        for i, j, value in operator.entries:
            contributions[j][classes[i]] += value
        signatures = [(classes[j], tuple(contributions[j])) for j in range(operator.n)]
        refined = partition(signatures)
        if refined == classes:
            break
        rounds.append({"before": groups, "after": max(refined)+1})
        classes = refined
    E = [[int(classes[j] == c) for j in range(operator.n)] for c in range(max(classes)+1)]
    return {"classes": list(classes), "E": E, "rounds": rounds,
            "claim": "closed aggregate partition; minimal among refinements of initial labels"}
