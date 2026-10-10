"""Optional NumPy benchmark: exact aggregate reduction of a periodic stencil.

Includes encoding and compilation amortization, and every-time output errors.
This is a structured demonstration, not an industrial performance claim.
"""
import argparse
import json
import os
import platform
import sys
import time
from pathlib import Path
from statistics import median
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from dynacomp.sparse import SparseOperator, compile_sparse, replay_sparse


def main():
    import numpy as np
    parser = argparse.ArgumentParser()
    parser.add_argument("--n", type=int, default=2048)
    parser.add_argument("--steps", type=int, default=200)
    parser.add_argument("--queries", type=int, default=256)
    args = parser.parse_args()
    if min(args.n, args.steps, args.queries) < 1:
        parser.error("Dimensions and budgets must be positive")
    n, steps, queries = args.n, args.steps, args.queries
    entries = [(i, i, "49/100") for i in range(n)]
    entries += [(i, (i+1)%n, "1/4") for i in range(n)]
    entries += [(i, (i-1)%n, "1/4") for i in range(n)]
    start = time.perf_counter()
    op = SparseOperator.from_entries(n, entries)
    packet = compile_sparse(op, [[0]]*n, [[1]*n])
    assert replay_sparse(packet)
    compilation = time.perf_counter()-start
    factor = float(packet["G"][0][0].split('/')[0])/float(packet["G"][0][0].split('/')[1])
    rng = np.random.default_rng(8128)
    initial = rng.normal(size=(queries, n))
    def full():
        x = initial.copy()
        y = [x.sum(axis=1)]
        for _ in range(steps):
            x = .49*x + .25*np.roll(x, 1, axis=1)+.25*np.roll(x, -1, axis=1)
            y.append(x.sum(axis=1))
        return np.asarray(y)
    def reduced():
        z = initial.sum(axis=1)
        y = [z.copy()]
        for _ in range(steps):
            z = factor*z
            y.append(z.copy())
        return np.asarray(y)
    # Warm both paths, then alternate measurements to reduce order bias.
    full(); reduced()
    timings = {"full": [], "reduced": []}
    outputs = {}
    for trial in range(3):
        methods = (("full", full), ("reduced", reduced))
        if trial%2:
            methods = methods[::-1]
        for name, function in methods:
            start = time.perf_counter()
            outputs[name] = function()
            timings[name].append(time.perf_counter()-start)
    f, r = median(timings["full"]), median(timings["reduced"])
    result = {"benchmark": "periodic stencil aggregate", "n": n, "rank": len(packet["E"]),
              "steps": steps, "queries": queries, "compile_and_replay_seconds": compilation,
              "full_batch_seconds": f, "reduced_batch_seconds_including_encoding": r,
              "execution_speedup": f/r, "first_batch_speedup_including_compilation": f/(r+compilation),
              "break_even_queries_estimate": compilation/((f-r)/queries) if f > r else None,
              "max_absolute_output_error_all_times": float(np.max(np.abs(outputs["full"]-outputs["reduced"]))),
              "raw_timings": timings, "python": platform.python_version(), "numpy": np.__version__,
              "platform": platform.platform(), "cpu_count": os.cpu_count(),
              "scope": "Exactly closed aggregate; not full-state reconstruction or generic simulation speedup"}
    (ROOT / "receipts" / "ring_benchmark.json").write_text(json.dumps(result, indent=2)+"\n")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
