# Changelog

## 0.2.0 — 10 October 2026

- `realization`: minimal input-output realization for zero initial state, with a Hankel-block minimality certificate checked against the original Markov parameters.
- `switched`: closure under every mode sequence, with experiment-labelled encoder rows; weighted automaton and HMM minimization (Schützenberger) with a Hankel certificate; equivalence checking that returns a shortest distinguishing word or an invariant-span certificate.
- `parametric` and `upoly`: generic observable, reachable and minimal orders for one real parameter, plus exact rank-drop loci (Cauchy–Binet, Faddeev–LeVerrier, Sturm isolation, rational roots found without factoring).
- `enrich`: certified greedy enrichment of approximate linear reductions to a target tolerance on supplied scenarios.
- `polynomial`: controlled polynomial closure, producing exact lifts z' = Σ uᵐ G_m z.
- `krylov`: a labelled breadth-first Krylov span with an incremental echelon form. The sparse compiler now uses it, roughly 28× faster on a 256-state stencil, and the 0.1 receipts regenerate byte for byte.
- New CLI kinds: `controlled_polynomial`, `realization`, `switched`, `automaton`, `equivalence`, `parametric` and `enrichment`. There are eight new example specifications and a heated-plate application.
- The suite grows from 38 to 60 tests.

## 0.1.0 — 10 October 2026

Initial implementation: rational linear and sparse observable closure, polynomial observable closure, finite partial-machine quotients, finite-horizon error envelopes, diagonal metric search, aggregate discovery, POD proposals and affine fast-forwarding.
