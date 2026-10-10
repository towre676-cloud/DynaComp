# DynaComp

**Compile the dynamics needed by a target, and account for what a reduction loses.** DynaComp constructs exact observable quotients of rational linear systems, sparse interaction networks, autonomous and controlled polynomial systems with finite observable closure, switched systems, and finite deterministic partial machines. It finds Hankel-certified minimal realizations, minimizes weighted automata and hidden Markov models, decides their equivalence with a shortest counterexample, and locates the exact parameter values where a reduction closes only by accident. For approximate linear models it computes rational finite-horizon output-error envelopes, and it can grow a reduced model until those envelopes meet a tolerance. This is a research implementation with executable demonstrations, not a general nonlinear reduction solver.

| Capability | Entry point | Accepted by |
|---|---|---|
| Observable quotient (dense / sparse) | `compile_linear`, `compile_sparse` | EA = GE, EB = H, DE = C |
| Minimal input-output realization | `minimal_realization` | AP = PF, B = PK, EF = GE, EK = H, DE = CP **and** a nonsingular Hankel block of the original |
| Switched-system closure | `compile_switched` | E A_s = G_s E for every mode s |
| Weighted automaton / HMM minimization | `minimize_automaton` | forward and backward identities plus a nonsingular Hankel block of word weights |
| Weighted language equivalence | `check_equivalence` | invariant span annihilating η₁ − η₂, or a shortest distinguishing word |
| Structural vs accidental closure, one parameter | `parametric_analysis` | Sturm-isolated real roots of the exact rank-drop polynomial |
| Certified greedy enrichment | `greedy_enrichment` | recomputed rational envelopes ≤ tolerance |
| Polynomial closure, autonomous / controlled | `compile_polynomial`, `compile_controlled_polynomial` | coefficient-by-coefficient pullback identities |
| Finite partial machines | `compile_finite` | homomorphism plus a distinguishing word per class pair |

The key question is whether a smaller state retains enough information to predict the requested output under future updates. When the initial representation fails, the linear compiler produces two states with the same representation and different next encoded states, then adds missing directions. Finite machines carry distinguishing experiments for every pair of retained classes. Polynomial compilation searches for features whose compositions with the update close exactly.

## Install and run

The exact core uses Python 3.10 or newer and the standard library. From the repository directory, install with `python -m pip install -e .`. NumPy is optional for snapshot proposals and the performance demonstration; install those extras with `python -m pip install -e '.[numerical]'`.

Run `python -m unittest discover -s tests -v` and `python examples/run_examples.py`. Compile a supplied model with `python -m dynacomp.cli compile examples/hidden_driver.json --output receipts/my_reduction.json`, then check its identities with `python -m dynacomp.cli replay receipts/my_reduction.json`. The CLI accepts `linear`, `sparse`, `polynomial`, `controlled_polynomial`, `finite`, `realization`, `switched`, `automaton`, `equivalence`, `parametric` and `enrichment` model kinds; every file in `examples/*.json` is a working specification. For example, `python -m dynacomp.cli compile examples/hmm_perturbed.json` returns the shortest observation word on which two weather HMMs disagree. Use rational strings such as `"1/2"`; the exact core rejects binary floating-point inputs.

## What works now

For linear systems, DynaComp constructs the smallest invariant row space containing the target rows. Its packet contains an encoder E, reduced transition G, reduced input H and output D satisfying **EA = GE, EB = H and DE = C**. These identities imply identical outputs for every initial state and every finite input sequence. Independent packet replay checks the identities without rerunning the discovery algorithm. The sparse backend computes row pullbacks from nonzero entries instead of constructing a dense transition matrix.

The polynomial backend discovers exact feature dynamics by repeated polynomial pullback and rational row-space reduction. The example `x_next = x + y², y_next = y` closes with the two features x and y². Degree, feature and iteration limits return a resource-limit result; they do not establish that a finite closure is impossible. The finite backend respects partial-action definedness and checks a separate distinguishing word for each pair of quotient classes.

Approximate models are handled through the residual identity for a rational lift V. Bounds include hidden initial-state error, transition residual and input residual, with target-aware propagation through powers of A. They remain valid on finite horizons even when A has transient growth or unstable modes. Positive diagonal metrics can supply an exact weighted infinity-norm contraction margin; failure to find such a metric does not prove instability. Optional POD snapshot proposals are candidates and require subsequent error analysis.

Sparse aggregate refinement discovers which groups have closed sum dynamics. A 16-site ring with one local defect refines into nine aggregate classes, preserving the reflection pairs. Exact affine block composition can fast-forward repeated input words in logarithmically many matrix compositions, subject to the cost of growing rational numerators and denominators.

## New in 0.2

Read [the 0.2 companion chapter](docs/EXTENSIONS.md) for derivations and limits.

- **Minimal realization with a minimality proof.** Reachable-then-observable quotient for zero initial state. Replay recomputes a Hankel block from the *original* Markov parameters, which shows that no smaller realization exists.
- **Switched systems, weighted automata and HMMs.** Closure under every mode sequence, with each encoder row labelled by its experiment word. Schützenberger minimization, and Tzeng-style equivalence that returns a *shortest* distinguishing word or an invariant-span certificate. Split-state weather HMMs are certified equivalent; a 1/50 emission perturbation is caught on a one-symbol word.
- **Parametric closure loci.** Generic observable, reachable and minimal orders over one real parameter, and *exactly* the parameter values where they drop, computed by Cauchy–Binet, Faddeev–LeVerrier and Sturm sequences without sampling. A mirror-symmetric rod is flagged as accidental exactly at t = 1/2.
- **Certified greedy enrichment.** Grows a lift from the target directions, choosing target-visible error snapshots, until the exact rational envelope of every scenario meets the tolerance. A six-cell chain is certified to 1/20 with five states.
- **Controlled polynomial lifts.** For x' = F(x, u) with arbitrary inputs, produces exact z' = Σ uᵐ G_m z.
- **Heated-plate application.** In a 12 × 12 plate whose target is the corner sum, the observable and Hankel-certified minimal orders are 17 of 144 states. Reduced evolution matches the full model exactly in rational arithmetic.
- **About 28× faster sparse compilation**, using an incremental echelon form with pivot-column coordinates.

## Demonstrated results

The suite now has **60 test methods**. They include brute-force word enumeration against automaton minimization and equivalence (including the shortest-counterexample claim), random Markov-parameter and Hankel-rank comparisons for realizations, pointwise rank checks against predicted parameter exceptions, and CLI round trips of every example. The initial 0.1 suite had 38 test methods, including the exhaustive 324 two-state partial-machine corpus, 60 random rational forced-system reductions, 30 random approximate-bound comparisons, packet-mutation rejection and trillion-step exact affine composition. A 64-site ring aggregate compiles to one state with multiplier 99/100. A stable nonnormal example exposes a first-step hidden-state output error of 8 and therefore rejects a tolerance-1 usable horizon beyond the initial time.

The optional benchmark in `examples/benchmark_ring.py` compiles a 2,048-state periodic stencil to its one-dimensional aggregate. For 256 queries of 200 steps, the retained local receipt records full execution at about 0.478 s, reduced execution including encoding at about 0.000422 s, and construction plus identity replay at about 0.121 s. The first batch is about **3.94× faster including construction**; the largest output difference over all sampled times is about 5.68e−14. This is an exactly closed aggregate demonstration. It does not reconstruct the full field or establish a generic simulation speedup. Rerun to obtain timings for your machine.

## Mathematical status and next development

Read [the monograph](docs/MONOGRAPH.md) for derivations, limitations and the redesign of the older IMP idea. [The 0.2 companion chapter](docs/EXTENSIONS.md) covers realization, switching, automata, parametric loci, enrichment and controlled lifts. [Development scope](docs/DEVELOPMENT.md) distinguishes implemented mechanisms from the planned nonlinear-memory, multiparameter and metric machinery. [Source notes](docs/SOURCES.md) identify the established methods this work builds on. Rational identity replay uses Python arithmetic; no Lean theorem or proof-kernel acceptance is claimed in this release.
