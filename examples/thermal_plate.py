"""Heated plate with corner sensors: exact sparse compilation and exact replay.

A 12 x 12 explicit heat stencil with uniform conductance 1/8 and ambient loss
1/100 is driven by a 2 x 2 heater in the center. The target is the summed
temperature of the four corner cells. The script compiles the observable
quotient from the sparse update, discovers the closed aggregate partition,
computes the input-output minimal order with a Hankel certificate, and then
checks reduced against full evolution in exact rational arithmetic.
"""
import json
import sys
import time
from fractions import Fraction
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from dynacomp import LinearModel
from dynacomp.sparse import SparseOperator, compile_sparse, replay_sparse
from dynacomp.structure import aggregate_partition
from dynacomp.realization import minimal_realization, replay_realization
from dynacomp.rational import matrix, mv, addv


def plate(N, conductance="1/8", loss="1/100"):
    k, entries = Fraction(conductance), []
    for r in range(N):
        for c in range(N):
            near = [(r+dr, c+dc) for dr, dc in ((1, 0), (-1, 0), (0, 1), (0, -1))
                    if 0 <= r+dr < N and 0 <= c+dc < N]
            entries.append((r*N+c, r*N+c, 1 - k*len(near) - Fraction(loss)))
            entries += [(r*N+c, rr*N+cc, k) for rr, cc in near]
    return SparseOperator.from_entries(N*N, entries)


def main(N=12, steps=40, minimal=True):
    op, n = plate(N), N*N
    center, corners = {N//2-1, N//2}, {0, N-1}
    B = [[int(i//N in center and i % N in center)] for i in range(n)]
    C = [[int(i//N in corners and i % N in corners) for i in range(n)]]
    t = time.perf_counter()
    packet = compile_sparse(op, B, C, max_rank=n)
    assert replay_sparse(packet)
    compile_seconds = time.perf_counter()-t
    labels = [(B[i][0], C[0][i]) for i in range(n)]
    partition = aggregate_partition(op, labels)
    receipt = {"grid": f"{N}x{N}", "full_dimension": n,
               "observable_order": len(packet["E"]), "observable_ranks_by_depth": packet["ranks"],
               "aggregate_classes": len(partition["E"]),
               "sparse_compile_and_replay_seconds": round(compile_seconds, 3)}
    if minimal:
        dense = [[0]*n for _ in range(n)]
        for i, j, v in op.entries:
            dense[i][j] = v
        t = time.perf_counter()
        realization = minimal_realization(LinearModel.from_rows(dense, B, C))
        assert replay_realization(realization)
        receipt["minimal_order"] = realization["minimal_order"]
        receipt["minimal_realization_with_hankel_certificate_seconds"] = round(time.perf_counter()-t, 3)
    # Exact replay: arbitrary initial field, heater pulses; every output must agree exactly.
    E, G, H, D = (matrix(packet[k]) for k in ("E", "G", "H", "D"))
    x = tuple(Fraction((7*i) % 11, 10) for i in range(n))
    z = mv(E, x)
    pulses = [[Fraction(1 if (s//5) % 2 == 0 else 0)] for s in range(steps)]
    for s, u in enumerate(pulses):
        x = tuple(a+b for a, b in zip(op.step(x), mv(matrix(B), u)))
        z = addv(mv(G, z), mv(H, u))
        assert mv(D, z) == mv(matrix(C), x), f"divergence at step {s}"
    receipt["exact_agreement_steps"] = steps
    receipt["final_corner_sum"] = str(mv(D, z)[0])
    receipt["final_corner_sum_approximately"] = f"{float(mv(D, z)[0]):.12g}"
    receipt["claim"] = ("exact for every initial field and input sequence of this model; "
                        "the reduction preserves the corner sum only, not the field")
    return receipt


if __name__ == "__main__":
    result = main()
    (ROOT / "receipts" / "thermal_plate.json").write_text(json.dumps(result, indent=2)+"\n")
    print(json.dumps(result, indent=2))
