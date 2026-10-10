# Run a model through DynaComp

Install with `python -m pip install -e '.[numerical]'`. Native execution additionally needs a GCC-compatible C compiler; the retained release was tested on Linux with GCC. No compiler is needed to compile or replay the mathematical JSON certificates. The native binary is built locally from packaged C source with separate multiplication/addition, contraction disabled, and no fast-math.

Try the small complete workflow:

```sh
dynacomp-pilot evaluate applications/customer_pilot/two_state.json --output pilot --queries 256 --trials 7 --batches 1000
```

Open `pilot/report.html` to see the bounds, measured execution, setup costs, failures and projected amortization. Change the batch count to see when setup pays back. `--batches` projects costs from measured medians; it does not run that many batches. Raw timing samples remain in `results.json`. Retained dense external plans and workload specifications use `.json.gz`; the commands read those directly.

For a supplied model, use the `dynacomp.workload/1` JSON schema demonstrated in `two_state.json`. The model is discrete-time `x_next = Ax + Bu`, `y = Cx`, with an output at time zero and after every transition. Dense `A,B,C` or the existing sparse model representation is accepted. Use rational strings, not JSON decimal numbers, for exact coefficients and domains. State the initial center/radius, one input radius row per transition, optional matching input centers, an absolute output error tolerance, and the rank search range. These boxes and the time horizon define the guarantee.

MAT data can be imported automatically:

```sh
dynacomp-pilot import-mat source.mat --dt 1/10 --output workload.json
dynacomp-pilot evaluate workload.json --output my-plan
```

Supplying `--dt` means the source is continuous-time and will be sampled with an augmented matrix exponential. Without it, matrices are treated as already discrete-time. Only identity mass matrix `E` and zero feedthrough `D` are supported. `--normalize-output` optionally divides all output rows by the largest row absolute sum, and records that divisor; default output units are preserved. Defaults are exploratory independent coordinate boxes and 20 steps: edit them to your actual intended operating domain. Imported guarantees start at the stored sampled coefficients; neither numerical sampling error nor physical calibration error is bounded.

Save the plan once, then execute your own batches:

```sh
dynacomp-pilot compile workload.json --output my-plan
dynacomp-pilot simulate my-plan --method observable --initial initial.npy --inputs inputs.npy --output outputs.npy
```

`initial.npy` has shape `(states, scenarios)` or `(states,)`; one initial field may be shared by the batch. `inputs.npy` has shape `(steps, channels, scenarios)`. Outputs have shape `(steps+1, outputs, scenarios)`. Loading uses `allow_pickle=False`. Preparation independently replays the saved certificate and verifies its source and domain binding. Every execution checks stored input values against inward-rounded endpoints of the exact boxes, rejects nonfinite inputs and the wrong horizon, and checks the floating-point environment.

Evaluate an already saved plan without searching again:

```sh
dynacomp-pilot run my-plan --queries 256 --trials 7 --batches 1000
```

For the external evidence used in this release, run `python benchmarks/customer_pilot.py`. It automatically ingests the retained, hashed sampled sources, searches both methods, and writes three complete reports under `receipts/customer_pilot`. Acquisition and sampling provenance is retained from the earlier external benchmark. To select one source, use `--case heat-cont`, `--case pde` or `--case build`. The original sources are pinned to the public SLICOT benchmark repository, rather than fetched from a changing URL at run time.

An accepted smaller model is not necessarily faster. The building example in this release is slower than full BLAS and has no measured break-even. Balanced truncation fails the declared box target within its searched ranks on that source. The reports retain both facts. The comparison is against our transparent SciPy square-root balanced implementation, not the compiled SLICOT toolbox. Both reduced methods run the same native engine and must pass the same combined model-and-execution error ceiling.
