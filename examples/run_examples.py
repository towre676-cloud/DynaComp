"""Regenerate deterministic mathematical receipts with standard-library Python."""
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "examples"))
from dynacomp.cli import compile_spec, replay
from dynacomp import LinearModel, finite_horizon_bound, usable_horizon
from dynacomp.metrics import discover_metric
from dynacomp.sparse import SparseOperator, compile_sparse, replay_sparse
from dynacomp.structure import aggregate_partition
from dynacomp.fastforward import repeat_input_block
from dynacomp import compile_linear
from dynacomp.parametric import classify
import thermal_plate


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
    for name in ("excitation", "switched_thermostat", "redundant_hmm", "hmm_equivalence",
                 "hmm_perturbed", "symmetric_rod", "bilinear_shear", "enrichment_chain"):
        packet = compile_spec(json.loads((ROOT / "examples" / f"{name}.json").read_text()))
        assert replay(packet)
        (out / f"{name}.json").write_text(json.dumps(packet, indent=2)+"\n")
        if name == "excitation":
            summary[name] = {"full": 3, "reachable": packet["reachable_dimension"],
                             "minimal_order": packet["minimal_order"], "hankel_certified": True}
        elif name == "switched_thermostat":
            summary[name] = {"full": 3, "reduced": len(packet["E"]), "experiments": packet["experiments"]}
        elif name == "redundant_hmm":
            summary[name] = {"original_order": packet["original_order"], "minimal_order": packet["order"]}
        elif name.startswith("hmm_"):
            summary[name] = {"equivalent": packet["equivalent"], **({"shortest_word": packet["word"],
                             "weights": packet["weights"]} if not packet["equivalent"] else {})}
        elif name == "symmetric_rod":
            summary[name] = {"generic_order": packet["minimal"]["generic_rank"],
                             "accidental_values": packet["minimal"]["rational_exceptions"],
                             "t=1/2": classify(packet, "1/2", "minimal"), "t=1/3": classify(packet, "1/3", "minimal")}
        elif name == "bilinear_shear":
            summary[name] = {"status": packet["status"], "features": len(packet["features"]),
                             "input_blocks": [b["input_exponent"] for b in packet["blocks"]]}
        else:
            summary[name] = {"status": packet["status"], "full": packet["full_dimension"],
                             "reduced": packet["reduced_dimension"],
                             "worst_bound_by_dimension": {h["dimension"]: h["approximately"] for h in packet["history"]}}
    plate = thermal_plate.main()
    plate.pop("sparse_compile_and_replay_seconds")
    plate.pop("minimal_realization_with_hankel_certificate_seconds")
    (out / "thermal_plate.json").write_text(json.dumps(plate, indent=2)+"\n")
    summary["thermal_plate"] = {k: plate[k] for k in ("full_dimension", "observable_order",
                                                      "aggregate_classes", "minimal_order", "exact_agreement_steps")}
    for name in ("uniform_directions", "sparse_thermal_family", "robust_switching", "validated_box"):
        packet = compile_spec(json.loads((ROOT / "examples" / f"{name}.json").read_text()))
        assert replay(packet)
        (out / f"{name}.json").write_text(json.dumps(packet, indent=2)+"\n")
        summary[name] = {"accepted": True, "schema": packet["schema"], "rank": packet["rank"]}
    (out / "summary.json").write_text(json.dumps(summary, indent=2)+"\n")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
