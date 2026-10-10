# Development scope

The current implementation is a target-driven dynamics compiler with four exact model backends, rational finite-horizon linear error envelopes, diagonal metric search, sparse aggregate discovery, rationalized POD proposals and affine block fast-forwarding. The oldest IMP ambition is retained as a question about sufficient state, rather than as an assertion that an orthogonal singular-vector projection equals an iterative limit.

The next substantial step is parameter-dependent compilation. A reduced representation may close at one coefficient choice and fail nearby. Polynomial or rational parameter identities should separate structurally valid families from accidental specialization, and expose rank-change loci. Harsh-Lattice-style symmetry candidates belong here when they induce checkable identities in the actual update operator.

Adaptive approximate reduction should then use target-weighted residuals to choose new coordinates. At present exact linear row closure repairs representations, but supplied approximate lifts are not automatically enriched. A practical enrichment algorithm needs a cost model, held-out forcing scenarios, and a stopping criterion expressed in the target tolerance. Construction cost and online encoding must remain in every performance comparison.

Nonlinear memory is a separate open problem. A projection can hide variables whose future effects cannot be expressed as an autonomous function of the retained coordinates. Finite delays or recurrent auxiliary variables may repair some cases, but no universal finite memory length follows. The first implementation should use polynomial examples with derived closures or explicit failure witnesses before making claims about general nonlinear simulations.

The present diagonal metric search is intentionally limited. General rational Lyapunov matrices, sparse semidefinite methods, continuous-time generators and robustness against switching deserve their own mathematical acceptance routes. Numerical optimization can propose a metric; exact or validated numerical inequalities must establish its claimed margin.

An initial physical application should be a parameterized sparse thermal or lattice-wave model with specified outputs and forcing distributions. Comparisons should include unreduced sparse evolution and appropriate established reduced-order methods at matched error tolerance. The eventual graphene Hamiltonian application requires correct complex arithmetic, unitary dynamics and target-specific error analysis; those are not implemented here.

The release does not claim a new universal reduction theorem, an industrial speedup, a replacement for pyMOR, automatic discovery of physical laws, or a formal Lean proof. Its concrete contribution is a coherent executable starting point that connects exact closure, useful failure witnesses, structural aggregation and bounded approximate evolution.
