# DynaComp

**Compile the dynamics needed by a target, and account for what a reduction loses.** DynaComp constructs exact observable quotients of rational linear systems, sparse interaction networks, autonomous polynomial systems with finite observable closure, and finite deterministic partial machines. It also computes rational finite-horizon output-error envelopes for supplied approximate linear models. This is an initial research implementation with executable demonstrations, not a general nonlinear reduction solver.

The key question is whether a smaller state retains enough information to predict the requested output under future updates. When the initial representation fails, the linear compiler produces two states with the same representation and different next encoded states, then adds missing directions. Finite machines carry distinguishing experiments for every pair of retained classes. Polynomial compilation searches for features whose compositions with the update close exactly.

## Install and run

The exact core uses Python 3.10 or newer and the standard library. From the repository directory, install with `python -m pip install -e .`. NumPy is optional for snapshot proposals and the performance demonstration; install those extras with `python -m pip install -e '.[numerical]'`.

Run `python -m unittest discover -s tests -v` and `python examples/run_examples.py`. Compile a supplied model with `python -m dynacomp.cli compile examples/hidden_driver.json --output receipts/my_reduction.json`, then check its identities with `python -m dynacomp.cli replay receipts/my_reduction.json`. The CLI accepts `linear`, `sparse`, `polynomial`, and `finite` model kinds. Use rational strings such as `"1/2"`; the exact core rejects binary floating-point inputs.

## What works now

For linear systems, DynaComp constructs the smallest invariant row space containing the target rows. Its packet contains an encoder E, reduced transition G, reduced input H and output D satisfying **EA = GE, EB = H and DE = C**. These identities imply identical outputs for every initial state and every finite input sequence. Independent packet replay checks the identities without rerunning the discovery algorithm. The sparse backend computes row pullbacks from nonzero entries instead of constructing a dense transition matrix.

The polynomial backend discovers exact feature dynamics by repeated polynomial pullback and rational row-space reduction. The example `x_next = x + y², y_next = y` closes with the two features x and y². Degree, feature and iteration limits return a resource-limit result; they do not establish that a finite closure is impossible. The finite backend respects partial-action definedness and checks a separate distinguishing word for each pair of quotient classes.

Approximate models are handled through the residual identity for a rational lift V. Bounds include hidden initial-state error, transition residual and input residual, with target-aware propagation through powers of A. They remain valid on finite horizons even when A has transient growth or unstable modes. Positive diagonal metrics can supply an exact weighted infinity-norm contraction margin; failure to find such a metric does not prove instability. Optional POD snapshot proposals are candidates and require subsequent error analysis.

Sparse aggregate refinement discovers which groups have closed sum dynamics. A 16-site ring with one local defect refines into nine aggregate classes, preserving the reflection pairs. Exact affine block composition can fast-forward repeated input words in logarithmically many matrix compositions, subject to the cost of growing rational numerators and denominators.

## Demonstrated results

The initial suite has **38 test methods**, including the exhaustive 324 two-state partial-machine corpus, 60 random rational forced-system reductions, 30 random approximate-bound comparisons, packet-mutation rejection and trillion-step exact affine composition. A 64-site ring aggregate compiles to one state with multiplier 99/100. A stable nonnormal example exposes a first-step hidden-state output error of 8 and therefore rejects a tolerance-1 usable horizon beyond the initial time.

The optional benchmark in `examples/benchmark_ring.py` compiles a 2,048-state periodic stencil to its one-dimensional aggregate. For 256 queries of 200 steps, the retained local receipt records full execution at about 0.478 s, reduced execution including encoding at about 0.000422 s, and construction plus identity replay at about 0.121 s. The first batch is about **3.94× faster including construction**; the largest output difference over all sampled times is about 5.68e−14. This is an exactly closed aggregate demonstration. It does not reconstruct the full field or establish a generic simulation speedup. Rerun to obtain timings for your machine.

## Mathematical status and next development

Read [the monograph](docs/MONOGRAPH.md) for derivations, limitations and the redesign of the older IMP idea. [Development scope](docs/DEVELOPMENT.md) distinguishes implemented mechanisms from the planned adaptive nonlinear, memory and parameter-dependent machinery. [Source notes](docs/SOURCES.md) identify the established methods this work builds on. Rational identity replay uses Python arithmetic; no Lean theorem or proof-kernel acceptance is claimed in this release.
