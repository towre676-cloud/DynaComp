# Changelog

## 0.5.0 — 10 October 2026

Sharp box-uniform general-direction enrichment with exact output-error kernels and attaining counterexamples; sparse coefficient-family closure; sharp full-model output ranges; thermal scenario workbench and optional batched QR execution with explicit rounding limits. Scenario-bound propagation now caches target pullbacks instead of full matrix powers, preserving earlier receipts. New CLI kinds `uniform_enrichment` and `sparse_family`, and command `dynacomp-thermal`. The suite grows from 95 to 117 methods.

## 0.4.0 — 10 October 2026

- `realization`: minimal input-output realization for zero initial state, with a Hankel-block minimality certificate checked against the original Markov parameters.
- `switched`: closure under every mode sequence, with experiment-labelled encoder rows; weighted automaton and HMM minimization (Schützenberger) with a Hankel certificate; equivalence checking that returns a shortest distinguishing word or an invariant-span certificate.
- `parametric` and `upoly`: generic observable, reachable and minimal orders for one real parameter, plus exact rank-drop loci (Cauchy–Binet, Faddeev–LeVerrier, Sturm isolation, rational roots found without factoring).
- `enrich`: certified greedy enrichment of approximate linear reductions to a target tolerance on supplied scenarios.
- `polynomial`: controlled polynomial closure, producing exact lifts z' = Σ uᵐ G_m z.
- `krylov`: a labelled breadth-first Krylov span with an incremental echelon form. The sparse compiler now uses it, roughly 28× faster on a 256-state stencil, and the 0.1 receipts regenerate byte for byte.
- New CLI kinds: `controlled_polynomial`, `realization`, `switched`, `automaton`, `equivalence`, `parametric` and `enrichment`. There are eight new example specifications and a heated-plate application.
- A randomized test checks that switched closure agrees with `compile_family`.
- The suite grows from 72 to 95 tests.

## 0.3.0 — 10 October 2026

Inventory-network scenario application (`applications/network_scenarios`): compact aggregate plans, conservation accounting, optional vectorized batch execution and honest amortization benchmarks.

## 0.2.0 — 10 October 2026

Adaptive box-bounded coordinate enrichment, least-cost bound-accepted subset search, coefficient-family closure with affine parameter conditions, and exact scalar output recurrences with input memory.

## 0.1.0 — 10 October 2026

Initial implementation: rational linear and sparse observable closure, polynomial observable closure, finite partial-machine quotients, finite-horizon error envelopes, diagonal metric search, aggregate discovery, POD proposals and affine fast-forwarding.
