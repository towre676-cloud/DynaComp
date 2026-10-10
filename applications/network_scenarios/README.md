# Network scenario workbench

This application takes a fixed linear routing network and evaluates replenishment scenarios. It returns regional inventory, dispatch totals, minimum stock margins and stock-floor violation counts. The compiler discovers a grouping whose sum dynamics preserve every requested output for every initial state and every future replenishment sequence. It splits proposed groups when local routing or dispatch differences invalidate closure. This is a working application prototype, and its examples are synthetic.

## Run the included example

From the repository root, install with `python -m pip install -e .` and run `python applications/network_scenarios/run_demo.py`. The example contains 32 stock cells in four regions, heterogeneous internal routing, three replenishment scenarios and forty steps. It compiles to four aggregate states. The output plan and complete scenario report are saved under `receipts`.

For direct CLI use, the following commands create a reusable plan, a JSON report and a CSV summary. The standard backend uses only Python's standard library.

```bash
dynacomp-network compile applications/network_scenarios/demo_network.json --output plan.json
dynacomp-network run plan.json applications/network_scenarios/demo_scenarios.json --output report.json --csv summary.csv
dynacomp-network run plan.json applications/network_scenarios/demo_scenarios.json --exact --output exact_report.json
dynacomp-network replay plan.json
```

Install `python -m pip install -e ".[numerical]"` and add `--backend numpy` to `run` for vectorized execution. Scenarios with equal horizons are evaluated together; different horizons are grouped automatically and results retain their supplied order. Exact arithmetic requires the standard backend. Ordinary floating-point results are not interval-certified, and stock-floor decisions close to equality should be resolved using exact arithmetic. The optional benchmark uses the same `NumpyEvolution` executor as this backend, rather than a separately optimized demonstration path.

## Model format and physical meaning

`demo_network.json` is a complete input example. A model has schema `dynacomp.network-model/1`, a positive node count `n`, unique `region_names`, one integer region index per node in `node_regions`, and sparse `transitions`. Each transition is `[destination, source, rational_fraction]`. Duplicate transitions are summed. Every supplied fraction must be nonnegative and each source column must sum to at most one. Integer values and rational or decimal strings are accepted; binary floating-point model values are rejected by the exact core.

The fraction missing from a source column is interpreted as completed dispatch from that source node's region. This interpretation is a modeling assumption: if the missing mass represents spoilage or leakage, the output should be interpreted accordingly. The current schema does not separate these different exit categories. Every node belongs to exactly one region, and every region contains at least one node.

`inputs` is a nonempty array of objects with unique `name` and a `distribution` array of `[node, weight]` pairs. Each distribution must sum exactly to one and have nonnegative weights. A scenario supplies one amount per named channel at each time. Those amounts are split over the nodes according to the declared weights; there is no hidden normalization or additional supply.

Stock is measured immediately before dispatch. At time t the current stock is routed and dispatched; replenishment u[t] arrives afterward, forming the next stock x[t+1] = A x[t] + B u[t]. Dispatch at t is calculated from x[t], not from the arriving supply. A horizon of T inputs therefore produces T+1 stock observations and T dispatch observations. Summed initial stock plus all supplied replenishment equals final stock plus cumulative dispatch in exact arithmetic. The report includes the conservation residual; numerical execution can accumulate small rounding differences.

This is a fixed-fraction linear model. It does not impose storage capacities, backlog, lead-time distributions, demand-limited dispatch or congestion. Stock floors are diagnostics only and do not change routing. Modeling a real warehouse requires agreeing that this operating regime fits its actual behavior. Existing nonlinear simulations can provide linear submodels, but this release does not construct those approximations or justify them from operational data.

## Scenario format and outputs

`demo_scenarios.json` uses schema `dynacomp.network-scenarios/1` and a nonempty `scenarios` array. Each item supplies a unique `name`, `initial_stock` with one nonnegative quantity per original node, and `replenishment` containing a time-ordered row per step, with one nonnegative quantity per input channel. Optional `stock_floor` contains one nonnegative floor per region. Exact scenario quantities follow the same integer/string convention as the model. An empty replenishment sequence is allowed and reports initial stock without any dispatch steps.

JSON reports preserve scenario order and list the region/channel names. Each scenario contains its complete stock and dispatch traces, cumulative dispatch, final stock, minimum margin above each floor, the number of times stock falls strictly below that floor, and its mass-balance residual. CSV summaries contain one row per scenario and region. Source and scenario fingerprints identify the precise JSON inputs used; they are provenance identifiers, not authentication signatures.

## Why compilation is trustworthy within this model

Let E sum the discovered classes. All nodes in a class must belong to the same reporting region, have the same dispatch fraction, and contribute the same total transition fraction to each destination class. These rational column conditions imply EA = GE. Summing each supply distribution over classes gives EB = H. Class-region indicators and dispatch fractions give DE = C. Consequently z[0] = E x[0], z[t+1] = G z[t] + H u[t], and y[t] = D z[t] agree exactly with the full model for every future input word.

The compact plan stores the class assignment rather than a dense binary encoder, plus the original model and rational G, H and D. Replay recomputes all column conditions and input/output coefficients without rerunning partition discovery. It accepts exact closure, not an independent proof that the saved refinement trace is minimal. The discovery algorithm itself computes the coarsest aggregate refinement of the region/dispatch labels. A `--max-groups` budget stops discovery before producing an over-budget encoder; this is a resource limit, not evidence that every possible reduction is impossible.

## Reproduce the performance evidence

Run `python applications/network_scenarios/benchmark.py`. The default cases have 128, 2,048 and 8,192 nodes, each with 128 scenarios and 96 steps. A fourth case has unique per-cell dispatch rates and retains all 32 nodes. The full baseline uses SciPy compressed sparse multiplication, and the structured cases also have an independent, correctly hand-built four-compartment baseline. Every path returns stock and dispatch at every time. Measurements use five warmed trials with rotated execution order and retain every raw timing.

The benchmark separately records rational construction/replay, numerical preparation, per-batch execution including initial encoding, first-batch cost and estimated amortization. Synthetic model generation, JSON I/O, scenario ingestion and installation are excluded. Thus the execution speedup is a kernel comparison, not an end-to-end CLI latency claim. The first batch may lose even when repeated execution is much faster. A hand aggregate can be faster still: the compiler's benefit is automated discovery and exact checking, not an inherently faster recurrence.
