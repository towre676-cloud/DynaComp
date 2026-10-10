"""Does a twenty-step thermal guarantee survive forever? Exact answer and repair.

The 0.5 workbench certifies a nine-state model of the 144-state heated plate:
corner-mean error below 0.1 for twenty steps, every initial field in [18, 22]
and every heater sequence in [0, 2]. This script asks the all-time question
with the same boxes held forever:

1. the nine-state model is REFUTED just after its horizon, and its all-time
   worst error is bracketed tightly near 3.27 (an attained scenario exists);
2. all-time enrichment, which can add steady-state moment directions, finds a
   seven-state model certified below 0.1 for every future time.
"""
import json
import sys
import time
from fractions import Fraction
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from dynacomp.thermal import plate_model
from dynacomp.uniform import uniform_enrichment
from dynacomp.forever import all_time_bound, all_time_enrichment, replay_all_time
from dynacomp.zonotope import sharp_zonotope_bound, box_generators


def binary64_witness(model, cand, T, center=20, radius=2, a=1, b=1):
    """Independent plain-float check: build the sign-pattern witness at time T and simulate the full plate."""
    F = lambda m: [[float(Fraction(x)) for x in row] for row in m]
    _, E, G, H = (F(x) for x in cand)
    entries = [(i, j, float(v)) for i, j, v in model.operator.entries]
    n, r = model.operator.n, len(G)
    Bv, C = [float(row[0]) for row in model.B], [float(x) for x in model.C[0]]
    D = [sum(C[i]*float(Fraction(cand[0][i][k])) for i in range(n)) for k in range(r)]
    ca, dg, J = C[:], D[:], []
    for _ in range(T):
        J.append(sum(x*y for x, y in zip(ca, Bv)) - sum(dg[k]*H[k][0] for k in range(r)))
        nxt = [0.0]*n
        for i, j, v in entries:
            nxt[j] += ca[i]*v
        ca, dg = nxt, [sum(dg[k]*G[k][c] for k in range(r)) for c in range(r)]
    K = [x - sum(dg[k]*E[k][c] for k in range(r)) for c, x in enumerate(ca)]
    sign = 1 if sum(k*center for k in K) + sum(J)*a >= 0 else -1
    x = [center + sign*(radius if k > 0 else -radius) for k in K]
    z = [sum(E[k][c]*x[c] for c in range(n)) for k in range(r)]
    for j in range(T):
        u = a + sign*(b if J[T-1-j] > 0 else -b)
        y = [0.0]*n
        for i, jj, v in entries:
            y[i] += v*x[jj]
        x = [p + q*u for p, q in zip(y, Bv)]
        z = [sum(G[k][c]*z[c] for c in range(r)) + H[k][0]*u for k in range(r)]
    return abs(sum(p*q for p, q in zip(C, x)) - sum(p*q for p, q in zip(D, z)))


def main(write=True):
    model = plate_model(12)
    n = model.operator.n
    center, radius, a, b, tol = [20]*n, [2]*n, [1], [1], "1/10"
    timings = {}
    t0 = time.perf_counter()
    finite = uniform_enrichment(model, center, radius, [b]*20, tol, input_centers=[a]*20, strategy="prefix")
    cand = (finite["V"], finite["E"], finite["G"], finite["H"])
    refuted = all_time_bound(model, *cand, center, radius, a, b, tolerance=tol)
    timings["refutation_s"] = time.perf_counter() - t0
    t0 = time.perf_counter()
    bracket = all_time_bound(model, *cand, center, radius, a, b)
    timings["converged_bracket_s"] = time.perf_counter() - t0
    t0 = time.perf_counter()
    forever = all_time_enrichment(model, center, radius, a, b, tol)
    timings["all_time_enrichment_s"] = time.perf_counter() - t0
    t0 = time.perf_counter()
    assert replay_all_time(forever)
    timings["independent_replay_s"] = time.perf_counter() - t0
    float_check = binary64_witness(model, cand, bracket["worst_time_lower_bound"])
    assert abs(float_check - float(Fraction(bracket["lower"]))) < 1e-6
    # Correlated initial uncertainty for the finite-horizon model (sharp, attained).
    T = 20
    independent = sharp_zonotope_bound(model, *cand, [0]*n, [[0]]*T, box_generators(n, radius))
    shared = sharp_zonotope_bound(model, *cand, [0]*n, [[0]]*T, [{"x0": {k: 2 for k in range(n)}}])
    summary = {
        "schema": "dynacomp.thermal-all-time/1",
        "domain": "initial field in [18,22]^144, heater in [0,2] at every step forever, corner-mean target",
        "finite_horizon_model": {"rank": finite["rank"], "certified_steps": 20,
                                 "twenty_step_sharp_bound": f'{float(Fraction(finite["certificate"]["worst_bound"])):.9g}'},
        "all_time_verdict_for_finite_model": {
            "status": refuted["status"], "violation_time": refuted["refutation"]["time"],
            "error_at_least": f'{float(Fraction(refuted["refutation"]["error_at_least"])):.9g}',
            "steady_state_error_at_nominal_heater": refuted["steady_state_error_approximately"][0]},
        "all_time_bracket_for_finite_model": {
            "status": bracket["status"], "lower": bracket["lower_approximately"],
            "upper": bracket["upper_approximately"], "attained_near_time": bracket["worst_time_lower_bound"],
            "steps_examined": bracket["steps_examined"],
            "independent_binary64_witness_error": f"{float_check:.9g}"},
        "all_time_model": {"status": forever["status"], "rank": forever["reduced_order"],
                           "lower": forever["lower_approximately"], "upper": forever["upper_approximately"],
                           "history": forever["enrichment_history"]},
        "correlated_initial_uncertainty_20_steps": {
            "independent_cells_pm2": independent["worst_bound_approximately"],
            "shared_offset_pm2": shared["worst_bound_approximately"]},
        "claim": ("rigorous rational/dyadic brackets; certification covers every t >= 0 for the declared boxes; "
                  "the model is synthetic and nothing is claimed about a physical plate"),
    }
    if write:
        (ROOT / "receipts" / "thermal_all_time.json").write_text(json.dumps(summary, indent=2) + "\n")
        (ROOT / "receipts" / "thermal_all_time_certificate.json").write_text(json.dumps(forever, indent=2) + "\n")
    return summary, timings


if __name__ == "__main__":
    summary, timings = main()
    print(json.dumps(summary, indent=2))
    print(json.dumps({k: round(v, 2) for k, v in timings.items()}, indent=2))
