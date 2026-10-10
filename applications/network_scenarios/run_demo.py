"""Reproduce complete scenario reports without optional dependencies."""
import json
import sys
from pathlib import Path
ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT/"src"))
from dynacomp.network import compile_network, run_scenarios


def main():
    here = Path(__file__).resolve().parent
    model = json.loads((here/"demo_network.json").read_text())
    scenarios = json.loads((here/"demo_scenarios.json").read_text())
    plan = compile_network(model)
    reduced = run_scenarios(plan, scenarios)
    full = run_scenarios(plan, scenarios, "full")
    maximum = max(abs(a-b) for s,t in zip(reduced["scenarios"],full["scenarios"]) for row,other in zip(s["stock_trace"]+s["dispatch_trace"],t["stock_trace"]+t["dispatch_trace"]) for a,b in zip(row,other))
    assert maximum < 1e-9
    (ROOT/"receipts"/"network_demo_plan.json").write_text(json.dumps(plan, indent=2)+"\n")
    (ROOT/"receipts"/"network_demo_report.json").write_text(json.dumps(reduced, indent=2)+"\n")
    print(json.dumps({"full_states": model["n"], "compiled_states": plan["rank"], "scenarios": len(reduced["scenarios"]), "max_absolute_difference": maximum}))


if __name__ == "__main__":
    main()
