"""Regenerate deterministic mathematical receipts with standard-library Python."""
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from dynacomp.cli import compile_spec, replay
from dynacomp import LinearModel, finite_horizon_bound, usable_horizon
from dynacomp.metrics import discover_metric
from dynacomp.sparse import SparseOperator, compile_sparse, replay_sparse
from dynacomp.structure import aggregate_partition
from dynacomp.fastforward import repeat_input_block
from dynacomp import compile_linear


def main():
    out = ROOT / "receipts"
    out.mkdir(exist_ok=True)
    summary = {}
    for name in ("hidden_driver", "nonlinear_shear", "partial_machine"):
        packet = compile_spec(json.loads((ROOT / "examples" / f"{name}.json").read_text()))
        assert replay(packet)
        (out / f"{name}.json").write_text(json.dumps(packet, indent=2)+"\n")
        summary[name] = {"accepted": True, "schema": packet["schema"]}
    n = 64
    entries = [(i, i, "49/100") for i in range(n)]
    entries += [(i, (i+1)%n, "1/4") for i in range(n)]
    entries += [(i, (i-1)%n, "1/4") for i in range(n)]
    op = SparseOperator.from_entries(n, entries)
    packet = compile_sparse(op, [[0]]*n, [[1]*n])
    assert replay_sparse(packet)
    (out / "ring64.json").write_text(json.dumps(packet, indent=2)+"\n")
    summary["ring"] = {"full_dimension": n, "reduced_dimension": len(packet["E"]), "factor": packet["G"][0][0]}
    model = LinearModel.from_rows([["1/2", 8], [0, "1/2"]], [[0], [0]], [[1, 0]])
    bound = finite_horizon_bound(model, [[1], [0]], [["1/2"]], [[0]],
                                 {"full": [0, 1], "reduced": [0]}, [[0]]*12)
    bound["description"] = "Stable eigenvalues do not suppress transient hidden-state output"
    bound["prefix_horizon_at_tolerance_1"] = usable_horizon(bound["bound"], 1)
    bound["metric"] = discover_metric(model.A)
    (out / "nonnormal_bound.json").write_text(json.dumps(bound, indent=2)+"\n")
    summary["nonnormal"] = {"first_step_error_bound": bound["bound"][1],
                            "prefix_horizon": bound["prefix_horizon_at_tolerance_1"],
                            "contractive_weighted_metric": bound["metric"]["contractive"]}
    defect_n = 16
    entries = [(i, i, "1/2") for i in range(defect_n)]
    entries += [(i, (i+1)%defect_n, "1/4") for i in range(defect_n)]
    entries += [(i, (i-1)%defect_n, "1/4") for i in range(defect_n)]
    entries += [(0, 0, "1/100")]
    discovered = aggregate_partition(SparseOperator.from_entries(defect_n, entries), [0]*defect_n)
    (out / "defect_partition.json").write_text(json.dumps(discovered, indent=2)+"\n")
    summary["defect_partition"] = {"full_dimension": defect_n, "aggregate_dimension": len(discovered["E"])}
    integrator = LinearModel.from_rows([[1, 1], [0, 1]], [[0], [1]], [[1, 0]])
    trillion = repeat_input_block(compile_linear(integrator), [2, 3], [[1]], 10**12)
    summary["fast_forward"] = {"steps": 10**12, "exact_output": str(trillion[0]), "method": "binary affine composition"}
    (out / "summary.json").write_text(json.dumps(summary, indent=2)+"\n")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
