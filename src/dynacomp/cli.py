"""JSON model compilation and independent identity replay."""
import argparse
import json
from pathlib import Path
from .linear import LinearModel, compile_linear, replay_linear
from .finite import FiniteModel, compile_finite, replay_finite
from .polynomial import compile_polynomial, unpack, replay_polynomial
from .sparse import SparseOperator, compile_sparse, replay_sparse
from .polynomial import compile_controlled_polynomial, replay_controlled_polynomial
from .realization import minimal_realization, replay_realization
from .switched import (compile_switched, replay_switched, WeightedAutomaton, minimize_automaton,
                       replay_automaton, check_equivalence, replay_equivalence)
from .parametric import parametric_analysis, replay_parametric
from .enrich import greedy_enrichment, replay_enrichment
from .families import compile_family, replay_family
from .memory import compile_output_memory, replay_memory
from .adaptive import adaptive_coordinates, replay_adaptive


def compile_spec(spec):
    kind = spec["kind"]
    if kind == "linear":
        return compile_linear(LinearModel.from_rows(spec["A"], spec["B"], spec["C"])).packet()
    if kind == "memory":
        return compile_output_memory(LinearModel.from_rows(spec["A"], spec["B"], spec["C"]))
    if kind == "family":
        return compile_family([LinearModel.from_rows(item["A"], item["B"], spec["C"]) for item in spec["coefficients"]])
    if kind == "adaptive":
        model = LinearModel.from_rows(spec["A"], spec["B"], spec["C"])
        return adaptive_coordinates(model, spec["initial_center"], spec["initial_radius"],
                                    spec["input_radii"], spec["tolerance"], **spec.get("options", {}))
    if kind == "finite":
        return compile_finite(FiniteModel(tuple(spec["observations"]), tuple(tuple(r) for r in spec["transitions"])))
    if kind == "polynomial":
        n = spec["variables"]
        update = [unpack(p, n) for p in spec["update"]]
        targets = [unpack(p, n) for p in spec["targets"]]
        if len(update) != n:
            raise ValueError("Update count must match variables")
        return compile_polynomial(update, targets, **spec.get("limits", {}))
    if kind == "sparse":
        op = SparseOperator.from_entries(spec["n"], spec["entries"])
        return compile_sparse(op, spec["B"], spec["C"], spec.get("max_rank", 128))
    if kind == "controlled_polynomial":
        n, k = spec["variables"], spec["inputs"]
        update = [unpack(p, n+k) for p in spec["update"]]
        targets = [unpack(p, n) for p in spec["targets"]]
        if len(update) != n:
            raise ValueError("Update count must match variables")
        return compile_controlled_polynomial(update, targets, k, **spec.get("limits", {}))
    if kind == "realization":
        return minimal_realization(LinearModel.from_rows(spec["A"], spec["B"], spec["C"]))
    if kind == "switched":
        return compile_switched(spec["modes"], spec["C"], spec.get("inputs"))
    if kind == "automaton":
        return minimize_automaton(WeightedAutomaton.from_spec(spec))
    if kind == "equivalence":
        return check_equivalence(WeightedAutomaton.from_spec(spec["first"]),
                                 WeightedAutomaton.from_spec(spec["second"]))
    if kind == "parametric":
        return parametric_analysis(spec["A"], spec["B"], spec["C"])
    if kind == "enrichment":
        model = LinearModel.from_rows(spec["A"], spec["B"], spec["C"])
        scenarios = [(s["x0"], s["inputs"]) for s in spec["scenarios"]]
        return greedy_enrichment(model, scenarios, spec["tolerance"], spec.get("max_dimension"))
    raise ValueError("Unknown model kind")


def replay(packet):
    readers = {"dynacomp.linear/1": replay_linear, "dynacomp.finite/1": replay_finite,
               "dynacomp.polynomial/1": replay_polynomial, "dynacomp.sparse/1": replay_sparse,
               "dynacomp.controlled_polynomial/1": replay_controlled_polynomial,
               "dynacomp.realization/1": replay_realization, "dynacomp.switched/1": replay_switched,
               "dynacomp.automaton/1": replay_automaton, "dynacomp.equivalence/1": replay_equivalence,
               "dynacomp.parametric/1": replay_parametric, "dynacomp.enrichment/1": replay_enrichment,
               "dynacomp.family/1": replay_family, "dynacomp.memory/1": replay_memory,
               "dynacomp.adaptive/1": replay_adaptive}
    if packet.get("schema") not in readers:
        raise ValueError("Unknown or unsuccessful packet")
    return readers[packet["schema"]](packet)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("operation", choices=["compile", "replay"])
    parser.add_argument("input", type=Path)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    try:
        data = json.loads(args.input.read_text())
        result = compile_spec(data) if args.operation == "compile" else {"accepted": replay(data)}
    except (ValueError, KeyError, TypeError, OverflowError, OSError) as exc:
        parser.exit(2, f"dynacomp: {exc}\n")
    text = json.dumps(result, indent=2) + "\n"
    if args.output:
        args.output.write_text(text)
    else:
        print(text, end="")


if __name__ == "__main__":
    main()
