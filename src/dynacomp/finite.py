"""Minimal all-future observational quotient for finite partial machines."""
from dataclasses import dataclass
from collections import deque


@dataclass(frozen=True)
class FiniteModel:
    observations: tuple
    transitions: tuple

    def __post_init__(self):
        n = len(self.observations)
        if not n or len(self.transitions) != n:
            raise ValueError("Need observations and transition row per state")
        m = len(self.transitions[0])
        if not m or any(len(row) != m for row in self.transitions):
            raise ValueError("Need a nonempty fixed action alphabet")
        if any(t is not None and (type(t) is not int or not 0 <= t < n)
               for row in self.transitions for t in row):
            raise ValueError("Successors must be states or None")


def partition(labels):
    distinct, result = [], []
    for label in labels:
        if label not in distinct:
            distinct.append(label)
        result.append(distinct.index(label))
    return tuple(result)


def distinguish(model, s, t):
    """Shortest experiment; definedness itself is observable."""
    queue, seen = deque([(s, t, ())]), set()
    while queue:
        s, t, word = queue.popleft()
        if (s, t) in seen:
            continue
        seen.add((s, t))
        if model.observations[s] != model.observations[t]:
            return list(word)
        for action, (a, b) in enumerate(zip(model.transitions[s], model.transitions[t])):
            if (a is None) != (b is None):
                return list(word + (action,))
            if a is not None:
                queue.append((a, b, word + (action,)))
    return None


def compile_finite(model):
    classes = partition(model.observations)
    rounds = []
    while True:
        signatures = [(classes[s], tuple(None if t is None else classes[t] for t in row))
                      for s, row in enumerate(model.transitions)]
        refined = partition(signatures)
        if refined == classes:
            break
        witnesses = []
        for s in range(len(classes)):
            for t in range(s):
                if classes[s] == classes[t] and refined[s] != refined[t]:
                    witnesses.append({"states": [t, s], "word": distinguish(model, t, s)})
        rounds.append({"before": max(classes)+1, "after": max(refined)+1, "witnesses": witnesses})
        classes = refined
    representatives = [classes.index(c) for c in range(max(classes)+1)]
    observations = [model.observations[s] for s in representatives]
    transitions = [[None if t is None else classes[t] for t in model.transitions[s]]
                   for s in representatives]
    witnesses = [{"classes": [a, b], "word": distinguish(model, representatives[a], representatives[b])}
                 for a in range(len(representatives)) for b in range(a)]
    return {"schema": "dynacomp.finite/1", "observations": list(model.observations),
            "transitions": [list(row) for row in model.transitions], "classes": list(classes),
            "quotient_observations": observations, "quotient_transitions": transitions,
            "distinguishing_witnesses": witnesses, "rounds": rounds}


def outcome(model, state, word):
    """Observe initial and every subsequent state, including undefined action."""
    trace = [model.observations[state]]
    for action in word:
        if type(action) is not int or not 0 <= action < len(model.transitions[0]):
            raise ValueError("Invalid action")
        state = model.transitions[state][action]
        if state is None:
            return (trace, False)
        trace.append(model.observations[state])
    return (trace, True)


def replay_finite(packet):
    """Check homomorphism and pair-specific witnesses, without partition refinement."""
    if packet.get("schema") != "dynacomp.finite/1":
        raise ValueError("Unknown schema")
    source = FiniteModel(tuple(packet["observations"]), tuple(tuple(r) for r in packet["transitions"]))
    target = FiniteModel(tuple(packet["quotient_observations"]), tuple(tuple(r) for r in packet["quotient_transitions"]))
    classes = packet["classes"]
    n = len(target.observations)
    if len(source.transitions[0]) != len(target.transitions[0]):
        raise ValueError("Action alphabets differ")
    if len(classes) != len(source.observations) or any(type(c) is not int or not 0 <= c < n for c in classes):
        raise ValueError("Invalid class map")
    if set(classes) != set(range(n)):
        raise ValueError("Class map is not onto")
    for s, c in enumerate(classes):
        if source.observations[s] != target.observations[c]:
            raise ValueError("Observation not preserved")
        for action, t in enumerate(source.transitions[s]):
            expected = None if t is None else classes[t]
            if expected != target.transitions[c][action]:
                raise ValueError("Transition/definedness not preserved")
    pairs = set()
    for witness in packet["distinguishing_witnesses"]:
        a, b = witness["classes"]
        if type(a) is not int or type(b) is not int or not 0 <= b < a < n:
            raise ValueError("Invalid witness pair")
        if outcome(target, a, witness["word"]) == outcome(target, b, witness["word"]):
            raise ValueError("Witness does not distinguish")
        pairs.add((a, b))
    if pairs != {(a, b) for a in range(n) for b in range(a)}:
        raise ValueError("Missing class-pair witnesses")
    return True
