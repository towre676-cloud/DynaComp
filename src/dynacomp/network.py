"""Inventory-network scenario compilation with exact closed regional outputs.

State is pre-dispatch stock; inputs arrive after dispatch each step. All column
mass absent from the transition is counted as dispatch from its source region.
Rates are fixed linear fractions, without demand caps or congestion feedback.
"""
import hashlib
import json
import math
from dataclasses import dataclass
from .rational import q, matrix, encoded, mv, addv
from .sparse import SparseOperator
from .structure import aggregate_partition


def fingerprint(spec):
    return hashlib.sha256(json.dumps(spec, sort_keys=True, separators=(",", ":"), allow_nan=False).encode()).hexdigest()


def _names(values, label):
    if not isinstance(values, list) or not values or any(not isinstance(v, str) or not v.strip() for v in values) or len(set(values)) != len(values):
        raise ValueError(f"Unique nonempty {label} names required")
    return tuple(values)


def _nonnegative(values, size, label):
    if not isinstance(values, (list, tuple)) or len(values) != size:
        raise ValueError(f"Wrong {label} dimension")
    result = tuple(q(v) for v in values)
    if any(v < 0 for v in result):
        raise ValueError(f"Nonnegative {label} required")
    return result


@dataclass(frozen=True)
class Network:
    op: SparseOperator
    regions: tuple
    names: tuple
    channels: tuple
    B: tuple
    exits: tuple

    @classmethod
    def load(cls, spec):
        if spec.get("schema") != "dynacomp.network-model/1":
            raise ValueError("Unknown network model schema")
        op = SparseOperator.from_entries(spec["n"], spec["transitions"])
        # Reject individual negative entries, even if duplicate summation cancels.
        if any(q(v) < 0 for _, _, v in spec["transitions"]):
            raise ValueError("Network transitions must be nonnegative")
        names = _names(spec["region_names"], "region")
        regions = tuple(spec["node_regions"])
        if len(regions) != op.n or any(type(v) is not int or not 0 <= v < len(names) for v in regions) or set(regions) != set(range(len(names))):
            raise ValueError("Every node needs a region and every region needs a node")
        totals = [q(0)]*op.n
        for _, j, v in op.entries:
            totals[j] += v
        if any(v > 1 for v in totals):
            raise ValueError("Transition columns must sum to at most one")
        inputs = spec["inputs"]
        channels = _names([item["name"] for item in inputs], "input")
        b = [[q(0)]*len(channels) for _ in range(op.n)]
        for k, item in enumerate(inputs):
            for node, weight in item["distribution"]:
                if type(node) is not int or not 0 <= node < op.n or q(weight) < 0:
                    raise ValueError("Invalid input distribution")
                b[node][k] += q(weight)
            if sum(row[k] for row in b) != 1:
                raise ValueError("Each input distribution must sum exactly to one")
        return cls(op, regions, names, channels, matrix(b), tuple(1-v for v in totals))

    def output(self, x):
        stock, dispatch = [q(0)]*len(self.names), [q(0)]*len(self.names)
        for i, value in enumerate(x):
            stock[self.regions[i]] += value
            dispatch[self.regions[i]] += self.exits[i]*value
        return tuple(stock+dispatch)


def _closed_coefficients(net, classes):
    """Check closure by column sums, independently of refinement discovery."""
    if len(classes) != net.op.n or any(type(c) is not int or c < 0 for c in classes):
        raise ValueError("Invalid aggregate classes")
    rank = max(classes)+1
    if set(classes) != set(range(rank)):
        raise ValueError("Aggregate classes must be consecutive and nonempty")
    representatives = [classes.index(c) for c in range(rank)]
    column = [[q(0)]*rank for _ in classes]
    for i, j, v in net.op.entries:
        column[j][classes[i]] += v
    for j in range(net.op.n):
        representative = representatives[classes[j]]
        if net.regions[j] != net.regions[representative] or net.exits[j] != net.exits[representative] or column[j] != column[representative]:
            raise ValueError("Aggregate does not preserve transition, stock or dispatch outputs")
    g = tuple(tuple(column[representatives[c]][d] for c in range(rank)) for d in range(rank))
    h = [[q(0)]*len(net.channels) for _ in range(rank)]
    for j, row in enumerate(net.B):
        for k, value in enumerate(row):
            h[classes[j]][k] += value
    d = [[q(0)]*rank for _ in range(2*len(net.names))]
    for c, j in enumerate(representatives):
        d[net.regions[j]][c] = q(1)
        d[len(net.names)+net.regions[j]][c] = net.exits[j]
    return g, matrix(h), matrix(d)


def compile_network(spec, max_groups=128):
    net = Network.load(spec)
    discovered = aggregate_partition(net.op, [(net.regions[j], net.exits[j]) for j in range(net.op.n)], max_groups)
    classes = discovered["classes"]
    g, h, d = _closed_coefficients(net, classes)
    return {"schema": "dynacomp.network-plan/1", "source_sha256": fingerprint(spec),
            "model": spec, "classes": classes, "G": encoded(g), "H": encoded(h), "D": encoded(d),
            "rank": len(g), "refinement_rounds": discovered["rounds"],
            "claim": "exact regional stock and dispatch for all initial states and input words; coarsest aggregate refinement of region/exit labels"}


def replay_network(packet):
    if packet.get("schema") != "dynacomp.network-plan/1" or fingerprint(packet["model"]) != packet["source_sha256"]:
        raise ValueError("Unknown plan schema or changed source model")
    net = Network.load(packet["model"])
    if type(packet["rank"]) is not int or packet["rank"] < 1 or packet["rank"] > net.op.n or len(packet["classes"]) != net.op.n or any(type(c) is not int or not 0 <= c < packet["rank"] for c in packet["classes"]):
        raise ValueError("Invalid plan rank or classes")
    g, h, d = _closed_coefficients(net, packet["classes"])
    if len(g) != packet["rank"] or any(matrix(packet[key]) != value for key, value in (("G", g), ("H", h), ("D", d))):
        raise ValueError("Network plan identities fail")
    return True


def load_scenarios(spec, net):
    if spec.get("schema") != "dynacomp.network-scenarios/1" or not isinstance(spec.get("scenarios"), list) or not spec["scenarios"]:
        raise ValueError("Nonempty scenario batch required")
    _names([s["name"] for s in spec["scenarios"]], "scenario")
    result = []
    for s in spec["scenarios"]:
        initial = _nonnegative(s["initial_stock"], net.op.n, "initial stock")
        if not isinstance(s["replenishment"], list):
            raise ValueError("Replenishment must be a time-ordered list")
        inputs = tuple(_nonnegative(row, len(net.channels), "input") for row in s["replenishment"])
        floor = _nonnegative(s.get("stock_floor", [0]*len(net.names)), len(net.names), "stock floor")
        result.append((s["name"], initial, inputs, floor))
    return result


def _trace(net, packet, initial, inputs, reduced, exact):
    convert = q if exact else float
    if reduced:
        x = [convert(0)]*packet["rank"]
        for j, value in enumerate(initial):
            x[packet["classes"][j]] += convert(value)
        g, h, d = (tuple(tuple(convert(q(v)) for v in row) for row in packet[k]) for k in ("G", "H", "D"))
        def dot(rows, values):
            return [sum(a*b for a, b in zip(row, values)) for row in rows]
        def output(state):
            return dot(d, state)
        def step(state, u):
            return [a+b for a, b in zip(dot(g, state), dot(h, u))]
    else:
        x = [convert(v) for v in initial]
        entries = [(i, j, convert(v)) for i, j, v in net.op.entries]
        b = [[convert(v) for v in row] for row in net.B]
        exits = [convert(v) for v in net.exits]
        def output(state):
            values = [convert(0)]*(2*len(net.names))
            for j, v in enumerate(state):
                values[net.regions[j]] += v
                values[len(net.names)+net.regions[j]] += exits[j]*v
            return values
        def step(state, u):
            result = [sum(a*v for a, v in zip(row, u)) for row in b]
            for i, j, v in entries:
                result[i] += v*state[j]
            return result
    trace = [output(x)]
    for u in inputs:
        x = step(x, [convert(v) for v in u])
        trace.append(output(x))
    if not exact and any(not math.isfinite(v) for row in trace for v in row):
        raise OverflowError("Nonfinite numerical output; rescale quantities or use exact arithmetic")
    return trace


def run_scenarios(packet, spec, engine="reduced", exact=False, backend="standard"):
    if backend not in ("standard", "numpy") or (exact and backend != "standard"):
        raise ValueError("Backend must be standard or numpy; exact arithmetic needs standard")
    if engine not in ("reduced", "full"):
        raise ValueError("Engine must be reduced or full")
    replay_network(packet)
    net = Network.load(packet["model"])
    results = []
    convert = str if exact else float
    scenarios = load_scenarios(spec, net)
    numerical_traces = {}
    if backend == "numpy":
        from .network_numeric import NumpyEvolution
        import numpy as np
        evolution = NumpyEvolution(net, packet, engine)
        buckets = {}
        for index, (_, initial, inputs, _) in enumerate(scenarios):
            buckets.setdefault(len(inputs), []).append(index)
        for horizon, indices in buckets.items():
            initial = np.asarray([[float(v) for v in scenarios[i][1]] for i in indices]).T
            inputs = np.asarray([[[float(v) for v in row] for row in scenarios[i][2]] for i in indices]).reshape(len(indices), horizon, len(net.channels)).transpose(1, 2, 0)
            output = evolution.simulate(initial, inputs)
            for column, index in enumerate(indices):
                numerical_traces[index] = output[:, :, column].tolist()
    for index, (name, initial, inputs, floor) in enumerate(scenarios):
        trace = numerical_traces[index] if backend == "numpy" else _trace(net, packet, initial, inputs, engine == "reduced", exact)
        r = len(net.names)
        cumulative = [sum(row[r+k] for row in trace[:-1]) for k in range(r)]
        supplied = sum(initial)+sum(sum(u) for u in inputs)
        residual = supplied-sum(trace[-1][:r])-sum(cumulative)
        margins = [[row[k]-(floor[k] if exact else float(floor[k])) for k in range(r)] for row in trace]
        results.append({"name": name, "horizon": len(inputs),
                        "stock_trace": [[convert(v) for v in row[:r]] for row in trace],
                        "dispatch_trace": [[convert(v) for v in row[r:]] for row in trace[:-1]],
                        "cumulative_dispatch": [convert(v) for v in cumulative],
                        "final_stock": [convert(v) for v in trace[-1][:r]],
                        "minimum_stock_margin": [convert(min(row[k] for row in margins)) for k in range(r)],
                        "stock_floor_violations": [sum(row[k] < 0 for row in margins) for k in range(r)],
                        "mass_balance_residual": convert(residual)})
    return {"schema": "dynacomp.network-report/1", "source_sha256": packet["source_sha256"],
            "scenario_sha256": fingerprint(spec), "engine": engine, "backend": backend, "arithmetic": "rational" if exact else "binary64",
            "full_states": net.op.n, "compiled_states": packet["rank"], "region_names": list(net.names),
            "input_names": list(net.channels), "scenarios": results,
            "scope": "fixed linear routing; floor checks are diagnostics, not capacity or demand feedback; binary64 output is not interval certified"}
