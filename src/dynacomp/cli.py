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
from .stability import (quadratic_certificate, replay_quadratic, continuous_certificate, replay_continuous,
                        common_quadratic_certificate, replay_common_quadratic)
from .forever import all_time_enrichment, replay_all_time
from .zonotope import sharp_zonotope_bound, replay_zonotope
from .uniform import _project
from .families import compile_family, replay_family
from .memory import compile_output_memory, replay_memory
from .adaptive import adaptive_coordinates, replay_adaptive
from .uniform import SparseLinearModel, uniform_enrichment, replay_uniform
from .sparse_family import compile_sparse_family, replay_sparse_family
from .robust import switching_enrichment, replay_switching
from .validated import validated_box,replay_validated
from .spectral import compile_spectral, compose_spectral, replay_spectral, replay_spectral_composition
from .ramanujan import ramanujan_value, replay_ramanujan_value
from .spectral import compress_spectral_modes, replay_spectral_subset
from .generating import RationalGF, from_recurrence, from_state, binomial_tail, replay_generating
from .binary64 import power_certificate, replay_power
from .machine import (certify_execution, compose_execution, replay_execution,
                      replay_machine_composition)


def compile_spec(spec):
    kind = spec["kind"]
    if kind == 'rational_gf':
        return RationalGF(spec['numerator'], spec['denominator']).packet()
    if kind == 'recurrence_gf':
        return from_recurrence(spec['coefficients'], spec['initial']).packet()
    if kind == 'state_gf':
        return from_state(spec['A'], spec['initial'], spec['C']).packet()
    if kind == 'binomial_tail':
        return binomial_tail(spec['alpha'], spec['argument'], spec['count'])
    if kind == 'rational_gf_tail':
        return RationalGF(spec['numerator'], spec['denominator']).tail(spec['argument'], spec['radius'], spec['count'])
    if kind == 'binary64_power':
        return power_certificate(spec['word'], spec['degree'])
    if kind == 'binary64_execution':
        return certify_execution(*(spec[k] for k in ('E','G','H','D')), spec['initial_center'],
                                 spec['initial_radius'], spec['input_radii'], spec['tolerance'], spec.get('input_centers'))
    if kind == 'machine_composition':
        return compose_execution(compile_spec(spec['source']), spec['tolerance'])
    if kind == "spectral_subset":
        return compress_spectral_modes(compile_spec(spec["source"]), spec["keep"], spec["tolerance"])
    if kind == "ramanujan_value":
        return ramanujan_value(spec["t"], spec.get("bits", 80))
    if kind == "spectral":
        return compile_spectral(spec["family"], spec["modes"], spec["tolerance"], spec.get("bits", 80))
    if kind == "spectral_composition":
        source = compile_spec(spec["source"])
        return compose_spectral(source, *(spec[k] for k in ("E", "G", "H", "D")),
                                spec["horizon"], spec["tolerance"], **spec.get("options", {}))
    if kind == "validated_box":
        return validated_box(SparseLinearModel.from_spec(spec['model']),spec['E'],spec['G'],spec['H'],spec['D'],spec['initial_center'],spec['initial_radius'],spec['input_radii'],spec['tolerance'],input_centers=spec.get('input_centers'),**spec.get('options',{}))
    if kind == "switching_enrichment":
        return switching_enrichment(spec["modes"], spec["C"], spec["initial_center"], spec["initial_radius"], spec["input_radii"], spec["tolerance"], input_centers=spec.get("input_centers"), **spec.get("options", {}))
    if kind == "uniform_enrichment":
        return uniform_enrichment(SparseLinearModel.from_spec(spec["model"]), spec["initial_center"], spec["initial_radius"], spec["input_radii"], spec["tolerance"], input_centers=spec.get("input_centers"), **spec.get("options", {}))
    if kind == "sparse_family":
        return compile_sparse_family(spec["modes"], spec["C"], spec.get("max_rank",128))
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
    if kind == "stability":
        if spec.get("time", "discrete") == "continuous":
            return continuous_certificate(spec["A"])
        return quadratic_certificate(spec["A"])
    if kind == "common_quadratic":
        return common_quadratic_certificate(spec["modes"])
    if kind == "all_time":
        return all_time_enrichment(SparseLinearModel.from_spec(spec["model"]), spec["initial_center"],
                                   spec["initial_radius"], spec["input_center"], spec["input_radius"],
                                   spec["tolerance"], spec.get("max_rank"))
    if kind == "zonotope":
        model = SparseLinearModel.from_spec(spec["model"])
        v, e, g, h, _ = _project(model, spec["encoder"])
        return sharp_zonotope_bound(model, v, e, g, h, spec["x0_center"], spec["input_centers"], spec["generators"])
    raise ValueError("Unknown model kind")


def replay(packet):
    if packet.get('schema') in ('exact.rational-gf/1','exact.rational-gf-tail/1','exact.binomial-tail/1'):
        return replay_generating(packet)
    if packet.get('schema') == 'exact.binary64-power/1':
        return replay_power(packet)
    if packet.get('schema') == 'dynacomp.binary64-execution/1':
        return replay_execution(packet)
    if packet.get('schema') == 'dynacomp.machine-composition/1':
        return replay_machine_composition(packet)
    if packet.get("schema") == "dynacomp.spectral-subset/1":
        return replay_spectral_subset(packet)
    if packet.get("schema") == "dynacomp.ramanujan-value/1":
        return replay_ramanujan_value(packet)
    if packet.get("schema") == "dynacomp.spectral/1":
        return replay_spectral(packet)
    if packet.get("schema") == "dynacomp.spectral-composition/1":
        return replay_spectral_composition(packet)
    readers = {"dynacomp.validated-box/1":replay_validated,"dynacomp.switching-enrichment/1": replay_switching, "dynacomp.uniform-enrichment/1": replay_uniform,
               "dynacomp.sparse-family/1": replay_sparse_family, "dynacomp.linear/1": replay_linear, "dynacomp.finite/1": replay_finite,
               "dynacomp.polynomial/1": replay_polynomial, "dynacomp.sparse/1": replay_sparse,
               "dynacomp.controlled_polynomial/1": replay_controlled_polynomial,
               "dynacomp.realization/1": replay_realization, "dynacomp.switched/1": replay_switched,
               "dynacomp.automaton/1": replay_automaton, "dynacomp.equivalence/1": replay_equivalence,
               "dynacomp.parametric/1": replay_parametric, "dynacomp.enrichment/1": replay_enrichment,
               "dynacomp.family/1": replay_family, "dynacomp.memory/1": replay_memory,
               "dynacomp.adaptive/1": replay_adaptive,
               "dynacomp.quadratic-metric/1": replay_quadratic, "dynacomp.continuous-metric/1": replay_continuous,
               "dynacomp.common-quadratic/1": replay_common_quadratic, "dynacomp.all-time-bound/1": replay_all_time,
               "dynacomp.zonotope-bound/1": replay_zonotope}
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
