# Switching scenario service

This application evaluates thermal scenarios while conductance changes separately for each scenario and time. Its rational plan covers every independent initial-temperature/heater box choice and every conductance value in [1/16,1/8] at each step. A nine-state approximation of the default 144-state plate stays below 0.1 corner-mean error through twenty steps. The certificate represents 1,048,576 endpoint mode words using 231 count states and retains an attaining experiment.

From the repository root, install `python -m pip install -e .` and run `python applications/switching_scenarios/run_demo.py`. It writes `receipts/thermal_switching_plan.json` and a readable summary. For supplied models, `dynacomp compile examples/robust_switching.json --output plan.json` and `dynacomp replay plan.json` exercise the new `switching_enrichment` kind. The specification has named sparse `modes`, shared `C`, initial/input centers and radii, tolerance, and options `max_rank`, `max_symbolic_states` and optional `seed_rows`.

The sharp switching route compresses words only when both full and projected operator families commute exactly. It retains order-dependent input effects by dynamic programming; simply comparing constant endpoint runs is insufficient. Noncommuting families are enumerated under an explicit budget. Resource exhaustion raises without accepting a partial bound. Exact common closure accepts zero error directly and can cover long horizons without either search.

Install `python -m pip install -e ".[numerical]"` for reusable numerical execution. Prepare the plan once, then use the same object for successive batches; the public API below includes domain checks by default.

```python
import json
import numpy as np
from dynacomp.prepared import PreparedSwitching

plan = json.load(open("receipts/thermal_switching_plan.json"))
service = PreparedSwitching(plan)
initial = np.full(144, 20.0)              # one field shared by every scenario
inputs = np.ones((20, 1, 256))            # time, channel, scenario
mixing = np.full((20, 256), 0.5)          # weight of the second named endpoint
outputs = service.simulate(initial, inputs, mixing)
print(outputs.shape, service.cache_stats)
```

Initial fields can also have shape [nodes, scenarios]. A shared field is encoded once per batch. A shared interpolation word can have shape [steps]. The plan is privately copied and checked once; later mutation of the original JSON object does not alter the executor. Reduced preparation replays the certificate; the full reference engine validates source/domain without charging it for a reduction proof. Numeric QR and evolution rounding are outside the rational model-reduction certificate.

Domain membership uses precomputed inward-directed binary64 endpoints, which give exact rational membership for finite binary64 quantities. This avoids per-value Fraction conversion and retains the rejection of rounded decimal values outside a rational endpoint. Setting `check_domain=False` is for already-validated workloads and removes the input-membership claim; the retained benchmark leaves checks enabled.

Run `python applications/switching_scenarios/benchmark.py` for three matched production cases: twenty-step robust approximate, twenty-step common exact, and four-hundred-step common exact. The long exact case runs about seven times faster including domain checks on the retained machine and amortizes construction/preparation after roughly eleven 256-scenario batches. Short cases need many more repetitions, and every first compiled batch loses after construction. The approximate plan is smaller but much more expensive to certify than the exact family. The benchmark separates all those costs and does not establish results on customer or measured physical data.

Derivations are in `docs/SWITCHING.md` and the complete monograph. Input weights are convex endpoint interpolation; arbitrary negative coefficient combinations are outside the approximate guarantee. The thermal stencil has fixed topology and spatially uniform conductance at each step. Uncertain bond-by-bond conductance, nonlinear heat laws and continuum error remain outside this application's contract.
