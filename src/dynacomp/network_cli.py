"""Compile inventory networks and report batches of replenishment scenarios."""
import argparse
import csv
import json
from pathlib import Path
from .network import compile_network, replay_network, run_scenarios


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="operation", required=True)
    compile_cmd = sub.add_parser("compile")
    compile_cmd.add_argument("model", type=Path)
    compile_cmd.add_argument("--max-groups", type=int, default=128)
    compile_cmd.add_argument("--output", type=Path, required=True)
    replay_cmd = sub.add_parser("replay")
    replay_cmd.add_argument("plan", type=Path)
    run_cmd = sub.add_parser("run")
    run_cmd.add_argument("plan", type=Path)
    run_cmd.add_argument("scenarios", type=Path)
    run_cmd.add_argument("--engine", choices=("reduced", "full"), default="reduced")
    run_cmd.add_argument("--exact", action="store_true")
    run_cmd.add_argument("--backend", choices=("standard", "numpy"), default="standard")
    run_cmd.add_argument("--output", type=Path, required=True)
    run_cmd.add_argument("--csv", type=Path)
    args = parser.parse_args()
    try:
        if args.operation == "compile":
            result = compile_network(json.loads(args.model.read_text()), args.max_groups)
        else:
            plan = json.loads(args.plan.read_text())
            if args.operation == "replay":
                result = {"accepted": replay_network(plan)}
            else:
                result = run_scenarios(plan, json.loads(args.scenarios.read_text()), args.engine, args.exact, args.backend)
        if args.operation == "run" and args.csv:
            with args.csv.open("w", newline="", encoding="utf-8") as stream:
                writer = csv.writer(stream)
                writer.writerow(["scenario", "region", "final_stock", "cumulative_dispatch", "minimum_stock_margin", "stock_floor_violations"])
                for s in result["scenarios"]:
                    for k, name in enumerate(result["region_names"]):
                        writer.writerow([s["name"], name]+[s[key][k] for key in ("final_stock", "cumulative_dispatch", "minimum_stock_margin", "stock_floor_violations")])
        text = json.dumps(result, indent=2, allow_nan=False)+"\n"
        if hasattr(args, "output"):
            args.output.write_text(text)
        else:
            print(text, end="")
    except (ValueError, KeyError, TypeError, OverflowError, ZeroDivisionError, OSError, ImportError) as exc:
        parser.exit(2, f"dynacomp-network: {exc}\n")


if __name__ == "__main__":
    main()
