"""Synthetic routing fixture, explicitly not measured warehouse data."""
import json
import random
from fractions import Fraction as Q
from pathlib import Path


def model(cells=8, defect=False, heterogeneous=False):
    if cells < 2:
        raise ValueError("At least two cells per region required")
    names = ["North", "East", "South", "West"]
    entries = []
    for region in range(4):
        base = region*cells
        for cell in range(cells):
            source = base+cell
            # Heterogeneous internal movement, but common external fractions.
            internal = Q(1+(cell*7)%19, 100)
            dispatch = Q(2+region, 100)
            if heterogeneous:
                dispatch += Q(cell+1, 100*cells)
            elif defect and region == 0 and cell == 0:
                dispatch += Q(1, 100)
            transfer = Q(3+region, 100)
            entries.extend([[source, source, str(1-internal-transfer-dispatch)],
                            [base+(cell+1)%cells, source, str(internal)],
                            [((region+1)%4)*cells+(cell*3+1)%cells, source, str(transfer)]])
    return {"schema": "dynacomp.network-model/1", "n": 4*cells,
            "region_names": names, "node_regions": [i//cells for i in range(4*cells)],
            "transitions": entries, "inputs": [{"name": name+" replenishment", "distribution": [[region*cells+j, str(Q(1,cells))] for j in range(cells)]} for region, name in enumerate(names)],
            "provenance": "Synthetic fixed-fraction inventory routing; not customer or measured data"}


def scenarios(cells=8, horizon=40, count=3, seed=20261010):
    rng = random.Random(seed)
    names = ["baseline", "North supply interruption", "South replenishment surge"]
    result = []
    for k in range(count):
        initial = [rng.randrange(5, 16) for _ in range(4*cells)]
        inputs = []
        for t in range(horizon):
            row = [max(0, cells//3+rng.randrange(0,4)) for _ in range(4)]
            if k%3 == 1 and horizon//4 <= t < 3*horizon//4:
                row[0] = 0
            if k%3 == 2 and t >= horizon//2:
                row[2] *= 3
            inputs.append(row)
        result.append({"name": names[k] if count == 3 else f"scenario_{k:04d}", "initial_stock": initial,
                       "replenishment": inputs, "stock_floor": [4*cells]*4})
    return {"schema": "dynacomp.network-scenarios/1", "scenarios": result}


if __name__ == "__main__":
    root = Path(__file__).resolve().parent
    for filename, data in (("demo_network.json", model()), ("demo_scenarios.json", scenarios())):
        (root/filename).write_text(json.dumps(data, indent=2)+"\n")
