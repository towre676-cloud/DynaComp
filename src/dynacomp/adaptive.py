"""Target-budgeted coordinate enrichment using exact rational box envelopes."""
from .rational import matrix, q, mul, encoded, mv
from .bounds import box_horizon_bound, usable_horizon
import heapq


def coordinate_candidate(model, selected):
    selected = tuple(selected)
    n = len(model.A)
    if not selected or len(set(selected)) != len(selected) or any(type(i) is not int or not 0 <= i < n for i in selected):
        raise ValueError("Need distinct valid coordinate indices")
    e = matrix([[int(j == i) for j in range(n)] for i in selected])
    v = matrix([[int(j == i) for i in selected] for j in range(n)])
    return {"selected": list(selected), "E": encoded(e), "V": encoded(v),
            "G": encoded(mul(mul(e, model.A), v)), "H": encoded(mul(e, model.B)),
            "D": encoded(mul(model.C, v))}


def adaptive_coordinates(model, initial_center, initial_radius, input_radii, tolerance,
                         initial_coordinates=None, max_rank=None, coordinate_costs=None):
    """Greedy target-error improvement per added-coordinate cost.

    Stopping is uniform over supplied boxes. Greedy selection is not globally
    rank-optimal, and costs are user-supplied estimates rather than timings.
    """
    n = len(model.A)
    tol = q(tolerance)
    if tol < 0:
        raise ValueError("Tolerance must be nonnegative")
    max_rank = n if max_rank is None else max_rank
    if type(max_rank) is not int or not 1 <= max_rank <= n:
        raise ValueError("Rank budget must lie between 1 and n")
    costs = [q(v) for v in (coordinate_costs if coordinate_costs is not None else [1]*n)]
    if len(costs) != n or any(v <= 0 for v in costs):
        raise ValueError("Coordinate costs must be strictly positive")
    selected = sorted(set(initial_coordinates if initial_coordinates is not None else
                          [j for j in range(n) if any(row[j] for row in model.C)]))
    if not selected or len(selected) > max_rank:
        raise ValueError("Initial coordinates exceed budget or are empty")
    # Validate once even when initial-coordinate indices are malformed.
    coordinate_candidate(model, selected)
    trace = []

    def assess(indices):
        candidate = coordinate_candidate(model, indices)
        bound = box_horizon_bound(model, candidate["V"], candidate["E"], candidate["G"], candidate["H"],
                                  initial_center, initial_radius, input_radii)
        worst = max(q(v) for v in bound["bound"])
        return candidate, bound, worst

    candidate, bound, worst = assess(selected)
    while True:
        trace.append({"selected": list(selected), "rank": len(selected), "worst_bound": str(worst),
                      "usable_horizon": usable_horizon(bound["bound"], tol)})
        if worst <= tol:
            status = "tolerance_met"
            break
        if len(selected) == max_rank:
            status = "rank_budget_exhausted"
            break
        trials = []
        for j in range(n):
            if j not in selected:
                c, b, w = assess(sorted(selected+[j]))
                score = (worst-w)/costs[j]
                trials.append((score, -w, -costs[j], -j, c, b, w))
        best = max(trials, key=lambda row: row[:4])
        _, _, _, neg_j, candidate, bound, worst = best
        selected = sorted(selected+[-neg_j])
        trace[-1]["chosen_coordinate"] = -neg_j
        trace[-1]["improvement_per_cost"] = str(best[0])
        trace[-1]["trials"] = [{"coordinate": -row[3], "worst_bound": str(row[6]),
                                 "improvement_per_cost": str(row[0])} for row in trials]
    return {"schema": "dynacomp.adaptive/1", "status": status, "tolerance": str(tol),
            "candidate": candidate, "box_bound": bound, "trace": trace,
            "cost": str(sum(costs[i] for i in selected)),
            "model": {k: encoded(getattr(model,k)) for k in ("A","B","C")},
            "domain": {"initial_center": [str(q(v)) for v in initial_center],
                       "initial_radius": [str(q(v)) for v in initial_radius],
                       "input_radii": [[str(q(v)) for v in u] for u in input_radii]},
            "max_rank": max_rank, "coordinate_costs": [str(v) for v in costs],
            "scope": "greedy coordinate family; uniform stated box and input budgets"}


def replay_adaptive(packet):
    """Replay the chosen candidate and bound, without rerunning greedy search."""
    from .linear import LinearModel
    if packet.get("schema") != "dynacomp.adaptive/1":
        raise ValueError("Unknown schema")
    model = LinearModel.from_rows(**packet["model"])
    candidate = coordinate_candidate(model, packet["candidate"]["selected"])
    if candidate != packet["candidate"]:
        raise ValueError("Candidate does not match selected coordinates")
    domain = packet["domain"]
    bound = box_horizon_bound(model, candidate["V"],candidate["E"],candidate["G"],candidate["H"],
                              domain["initial_center"],domain["initial_radius"],domain["input_radii"])
    if bound != packet["box_bound"]:
        raise ValueError("Uniform box bound does not match")
    tolerance = q(packet["tolerance"])
    budget = packet["max_rank"]
    costs = [q(v) for v in packet["coordinate_costs"]]
    if tolerance < 0 or type(budget) is not int or not 1 <= budget <= len(model.A):
        raise ValueError("Invalid tolerance/rank budget")
    if len(costs) != len(model.A) or any(v <= 0 for v in costs):
        raise ValueError("Invalid coordinate costs")
    if q(packet["cost"]) != sum(costs[i] for i in candidate["selected"]):
        raise ValueError("Reported cost does not match")
    rank = len(candidate["selected"])
    success = max(q(v) for v in bound["bound"]) <= tolerance
    if rank > budget or packet["status"] != ("tolerance_met" if success else "rank_budget_exhausted"):
        raise ValueError("Reported status does not follow from tolerance")
    if not success and rank != budget:
        raise ValueError("Rank budget not exhausted")
    return True


def optimal_coordinates(model, initial_center, initial_radius, input_radii, tolerance,
                        required_coordinates=None, max_rank=None, coordinate_costs=None,
                        max_evaluations=10000):
    """Best-first exact least-cost search within a specified coordinate family.

    The first feasible popped subset has minimum additive cost among subsets
    meeting this particular conservative box-envelope criterion. No monotonic
    error assumption or simulated-trace extrapolation is used.
    """
    n=len(model.A)
    required=tuple(sorted(set(required_coordinates if required_coordinates is not None else
                             [j for j in range(n) if any(row[j] for row in model.C)])))
    coordinate_candidate(model,required)
    limit=n if max_rank is None else max_rank
    if type(limit) is not int or not len(required)<=limit<=n:
        raise ValueError("Invalid rank limit")
    if type(max_evaluations) is not int or max_evaluations<1:
        raise ValueError("Positive evaluation budget required")
    tol=q(tolerance)
    costs=[q(v) for v in (coordinate_costs if coordinate_costs is not None else [1]*n)]
    if tol<0 or len(costs)!=n or any(v<=0 for v in costs):
        raise ValueError("Invalid tolerance or coordinate costs")
    missing=[i for i in range(n) if i not in required]
    heap=[(sum(costs[i] for i in required),len(required),(),0)]
    evaluations=[]
    while heap and len(evaluations)<max_evaluations:
        cost,rank,extra,start=heapq.heappop(heap)
        selected=sorted(required+extra)
        candidate=coordinate_candidate(model,selected)
        bound=box_horizon_bound(model,candidate['V'],candidate['E'],candidate['G'],candidate['H'],
                                initial_center,initial_radius,input_radii)
        worst=max(q(v) for v in bound['bound'])
        evaluations.append({'selected':selected,'cost':str(cost),'worst_bound':str(worst)})
        if worst<=tol:
            return {'schema':'dynacomp.coordinate-search/1','status':'optimal_for_bound_criterion',
                    'candidate':candidate,'box_bound':bound,'cost':str(cost),'evaluations':evaluations,
                    'criterion':'least additive cost among bound-accepted coordinate sets containing required coordinates',
                    'required_coordinates':list(required),'tolerance':str(tol),
                    'frontier_lower_cost':str(heap[0][0]) if heap else None}
        if rank<limit:
            for j in range(start,len(missing)):
                index=missing[j]
                heapq.heappush(heap,(cost+costs[index],rank+1,extra+(index,),j+1))
    return {'schema':'dynacomp.coordinate-search/1',
            'status':'evaluation_budget_exhausted' if heap else 'no_bound_accepted_subset',
            'evaluations':evaluations,'tolerance':str(tol),'required_coordinates':list(required),
            'frontier_lower_cost':str(heap[0][0]) if heap else None,
            'scope':'failure of bounded search or envelope criterion is not impossibility for the true dynamics'}
