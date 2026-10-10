# DynaComp 0.8: stable proposals and inexpensive rigorous box brackets

The external SLICOT benchmark in 0.7.1 exposed two concrete problems. Exact rational powers made heat-model construction and independent replay cost about 73 seconds. The generic orthogonal discovery heuristic gave the building-motion model an eight-state worst error exceeding twelve thousand. Version 0.8 addresses both on the same retained sampled sources, without changing their coordinate boxes, twenty-step horizon or normalized tolerance of 0.01. It combines Lyapunov-coordinate discovery, exact reduced contraction checks, bounded-precision integer sensitivity propagation, and reusable execution of the matrices actually accepted.

The other session's work is already integrated: 0.7 supplies exact Schur–Cohn/Hurwitz decisions, general rational quadratic metrics, all-time brackets for persistent boxes, moment enrichment and correlated finite-horizon domains. The seven-state plate result below 0.1 forever remains in `docs/FOREVER.md` and the complete monograph. The new external seven-state heat model is a different model and has a twenty-step guarantee. Neither its finite-horizon acceptance nor its reduced stability proves that its error tolerance holds forever.

## A stability-aware discovery route

The optional `StableProposalContext` proposes a numerical full metric P through the discrete Lyapunov equation Aᵀ P A − P = −I, factors P = LᵀL and works in coordinates y = Lx. A projection which is orthogonal in these coordinates is generally oblique in the original physical coordinates. Its reduced coordinates need not be original temperatures or displacements. The discovery starts with the transformed target C L⁻¹ and adds missing directions most visible to the requested output error, with a pullback fallback. The initial-state error is still scored in the original coordinate box, not a changed ellipsoid.

For a basis Q in the transformed coordinates, the proposal has E = QᵀL, G = QᵀL A L⁻¹Q, H = QᵀL B and D = C L⁻¹Q. The full metric and basis use binary64 and are proposals only. Each resulting E, G, H and D is bound as the exact dyadic rational represented by its binary64 value. Acceptance then checks **I − GᵀG ≻ 0** by exact rational symmetric elimination. This proves homogeneous reduced stability: every nonzero reduced state loses Euclidean energy, and a contraction factor below one exists. This check uses the actual supplied G. It does not require trusting the numerical full metric, a numerical eigenvalue or a claimed projection identity.

The test is sufficient, not necessary. A stable matrix may fail this Euclidean contraction check and admit a different quadratic metric. Version 0.7 provides the more general route. Failure of the full numerical Lyapunov/Cholesky proposal is not an exact instability decision. Conditioning and finite precision can interfere with discovery. The packet therefore distinguishes a proved reduced contraction, a refuted accuracy tolerance, an undecided accuracy bracket and a tolerance met without a requested stability proof. Stability of the reduced homogeneous system does not itself establish all-time accuracy or controller stability on the full plant.

The context computes its numerical metric once and extends a nested basis across requested ranks. Asking for a smaller rank later takes a prefix. The external search starts at rank four and grows through at most sixteen until both accuracy and requested reduced contraction are accepted. This is heuristic search, not minimum-rank optimization. In particular a failed rank-r proposal does not prove that every rank-r model must fail.

## Fixed-model error without a left-inverse prerequisite

The fast acceptance route takes a full source x′ = Ax + Bu, y = Cx and any supplied reduced system z′ = Gz + Hu, ŷ = Dz with initialization z₀ = Ex₀. At time t its error kernels are

K_t = C Aᵗ − D Gᵗ E, and J_l = C Aˡ B − D Gˡ H.

Thus y_t − ŷ_t = K_t x₀ + Σ_{j<t} J_{t−1−j} u_j. This identity requires matching dimensions, but **does not require EV = I**, an invariant subspace or a lift V. That makes the route useful for numerically assembled encoders, quantized models and standard balanced-truncation proposals. It certifies their actual supplied matrices rather than replacing them with a different idealized projection.

For independent boxes x₀ = c ± r and u_j = a_j ± b_j, each component's exact supremum is |K_t c + Σ J a_j| + Σ |K_t| r + Σ |J| b_j. The earlier exact module evaluates that formula with rational powers. The new module encloses it using outward dyadic rows. Domain centers and radii remain exact rational values. Output times and components are bounded separately, and the packet retains the worst upper bound across them.

## Integer propagation with an explicit rounding budget

With grid scale S = 2^bits, each coefficient or propagated row component is enclosed by a midpoint integer X and a nonnegative radius integer R, meaning the true quantity lies in [(X−R)/S, (X+R)/S]. Rational-to-grid conversion is performed exactly. A coefficient exactly on the grid has radius zero; other nearest-rounded coefficients have radius one. No ordinary floating-point multiplication participates in bound acceptance.

For row midpoint X_i and radius R_i, and matrix midpoint M_ij and radius T_ij, form the exact integer dot Z_j = Σ_i X_i M_ij and spread U_j = Σ_i (|X_i| T_ij + R_i (|M_ij| + T_ij)). Choose the nearest integer X′_j to Z_j/S, then set

R′_j = ceil((U_j + |Z_j − S X′_j|)/S).

The spread accounts for row uncertainty, coefficient uncertainty and their product. The last term accounts for the actual integer-dot rounding distance. Induction proves that every propagated component encloses its exact rational counterpart. Kernel subtraction adds radii. Exact dyadic cancellations can retain zero radius; no unconditional rounding allowance is injected when the calculation is exact.

Each kernel enclosure produces an upper support bound. A retained corner is chosen from the midpoint signs and the midpoint nominal-error orientation. Its signed support has a rigorous lower bound obtained by subtracting each uncertainty radius. A positive lower bound therefore identifies an actual admissible scenario whose absolute error reaches or exceeds it. A zero lower bound is trivial and does not claim sharp attainment. The result is a bracket around the sharp worst error, not a claim to have computed its exact rational value.

An upper bound below tolerance certifies the complete box domain. A witness lower bound above tolerance refutes the tolerance. When neither inequality resolves the question, the engine can refine through a declared sequence of precisions; the default is 64, 96, 128 and 192 bits. A exhausted precision budget gives `undecided`, never acceptance from an incomplete calculation. The chosen precision and bracket replay independently. Componentwise intervals can widen severely in long or nonnormal problems, so adaptive precision is not a universal cure for interval dependency or long-horizon cost. Centers and radii can also grow in magnitude; bounded grid denominators are not a constant-bit-complexity claim.

`BoxKernelCache` retains full-source pullbacks and forcing kernels separately by precision. Candidate search reuses them while recomputing reduced kernels. The cache checks source and horizon identity, and it never substitutes a bound from a different model. Independent packet replay constructs a fresh source calculation. It checks the chosen precision's bracket, witness, rank, acceptance status and any requested reduced contraction; descriptive discovery history is not a mathematical premise.

## External results

The sources are the exact stored sampled models from 0.7.1. Both the upstream MAT digests and canonical sampled-source digests are pinned in `benchmarks/slicot_manifest.json`. This avoids quietly changing discretization, output normalization or floating matrix-exponential results between comparisons. The workloads remain generated controlled scenarios on published engineering matrices, not new physical measurements. Complete receipts, rank-search histories, twelve compressed packets and raw timings are in `receipts/external_validated`.


For **heat-cont**, the earlier eight-state exact construction and replay cost 73.041 seconds. Rediscovery, dyadic acceptance and reusable setup including independent replay cost 1.209 seconds, a 60.4-fold cost reduction on the retained machine. The newly bound execution matrices are not bit-identical to the earlier exact Echelon-coordinate model, but represent its proposed subspace numerically and are accepted or rejected on their own coefficients. Their worst upper bound is 0.0003080475127, their bracket width 4.01e-16, and their status tolerance_met.

The stability-aware search chooses **7 states** with status **tolerance_met**, a rigorous upper error **0.006988834282**, and exact reduced contraction accepted. Its repeated production batches run **43.77×** faster than full execution. Discovery/certification plus prepared setup cost 1.352 seconds, and amortize after about **104 batches** of 256 scenarios. The matched-tolerance balanced search chooses 9 states with status tolerance_met, upper bound 0.0005995394883, execution speedup 39.74× and estimated amortization 114 batches. Both numerical methods receive the same dyadic acceptance machinery.

For **pde**, the earlier eight-state exact construction and replay cost 8.541 seconds. Rediscovery, dyadic acceptance and reusable setup including independent replay cost 0.300 seconds, a 28.4-fold cost reduction on the retained machine. The newly bound execution matrices are not bit-identical to the earlier exact Echelon-coordinate model, but represent its proposed subspace numerically and are accepted or rejected on their own coefficients. Their worst upper bound is 6.500117043e-07, their bracket width 1.12e-16, and their status tolerance_met.

The stability-aware search chooses **4 states** with status **tolerance_met**, a rigorous upper error **0.002163652598**, and exact reduced contraction accepted. Its repeated production batches run **15.51×** faster than full execution. Discovery/certification plus prepared setup cost 0.171 seconds, and amortize after about **50 batches** of 256 scenarios. The matched-tolerance balanced search chooses 4 states with status tolerance_met, upper bound 0.003059127706, execution speedup 15.60× and estimated amortization 52 batches. Both numerical methods receive the same dyadic acceptance machinery.

For **build**, the earlier eight-state exact construction and replay cost 1.978 seconds. Rediscovery, dyadic acceptance and reusable setup including independent replay cost 0.100 seconds, a 19.9-fold cost reduction on the retained machine. The newly bound execution matrices are not bit-identical to the earlier exact Echelon-coordinate model, but represent its proposed subspace numerically and are accepted or rejected on their own coefficients. Their worst upper bound is 12371.6251, their bracket width 2.51e-07, and their status refuted.

The stability-aware search chooses **13 states** with status **tolerance_met**, a rigorous upper error **0.0009013591345**, and exact reduced contraction accepted. Its repeated production batches run **3.35×** faster than full execution. Discovery/certification plus prepared setup cost 0.325 seconds, and amortize after about **490 batches** of 256 scenarios. The matched-tolerance balanced search chooses 16 states with status refuted, upper bound 40.58206631, execution speedup 2.95× and estimated amortization 424 batches. Both numerical methods receive the same dyadic acceptance machinery.

The building case changes qualitatively. The generic eight-state proposal remains refuted around error 12,371.6. In Lyapunov coordinates the twelve-state proposal is still refuted at about 0.0161, while thirteen states certify error below 0.000902. The witness reproduces the bracket numerically. Balanced discovery through sixteen states is still refuted under this initial-state box objective, although its responses from rest remain useful. These results justify stability-aware discovery for this case, not general superiority over balanced truncation or a proof of optimal rank.

The ordinary methods remain competitive. On the PDE model both searches accept four states, have nearly identical warm performance and roughly fifty-batch payback. On heat the stability-aware search accepts seven states; balanced search accepts nine and gives a smaller error bound. The stability-aware heat result meets the declared tolerance with fewer states and has somewhat faster execution, while a supplied twelve-state balanced model is especially cheap to discover. The expensive old certification step is improved for all of them. First compiled batches remain slower than full simulation, and timing differences between close methods should not be overinterpreted.

## Execution and reproduction

`PreparedValidated` independently replays a packet once, copies it privately, prepares its actual E, G, H and D, and caches all numerical matrices. Full execution validates only the source and domain; no reduction proof is charged to that baseline. Every batch checks shapes, horizon, finite data and box membership. Rational box endpoints are directed inward to binary64 once, making membership exact for finite binary64 batch values. A shared initial field is encoded once per batch. An unaccepted reduction is rejected by default; the benchmark explicitly enables it solely to execute retained counterexamples. There is no per-batch QR or reconstruction of a different reduced operator.

The mathematical bracket excludes coefficient conversion and numerical evolution rounding. In the retained benchmarks every new matrix coefficient is an exact binary64 dyadic, so coefficient loading recovers its proposed value. Generic rational packets may incur additional conversion rounding. Numerical outputs are compared with the rigorous model-reduction bound using an explicitly separate small rounding allowance; that allowance is not a certified numerical integration bound.

From the repository root install `python -m pip install -e ".[numerical]" threadpoolctl`, then run `OPENBLAS_NUM_THREADS=1 python benchmarks/external_validated.py`. On PowerShell set `$env:OPENBLAS_NUM_THREADS='1'` and run the Python command. The release includes the retained sampled-source packets, so this benchmark requires no new download. Rebuilding samples from upstream raw data on another numerical platform can change their exact digests; that is a new source experiment, not the identical pinned comparison. `benchmarks/README.md` explains upstream acquisition separately.

The timing uses one BLAS thread, five warmed alternating trials and identical output histories for every engine. Each workload contains random box members, persistent forcing from rest, a coherent initial field and each candidate's own worst corner. Source input/parsing is common and reported separately. Search construction and prepared setup including independent replay are included in amortization. Gzip serialization, JSON transfer and optional-library installation add further startup costs and are excluded. The first orthogonal candidate pays the first shared source-kernel pass. Both stability and balanced searches prepare their full metric or Gramian factorization once, rather than repeatedly charging full solves at every tested rank.

A small complete API route is:

```python
from dynacomp import SparseLinearModel, StableProposalContext, validated_box
from dynacomp.validated_numeric import PreparedValidated

source = SparseLinearModel.from_spec({
    'n': 2, 'entries': [[0, 0, '1/2'], [0, 1, 8], [1, 1, '1/2']],
    'B': [[0], [1]], 'C': [[1, 0]],
})
proposal = StableProposalContext(source, 20).propose(2)
packet = validated_box(source, *(proposal[k] for k in ('E','G','H','D')),
                       [0,0], [1,1], [[1]]*20, '1/1000000',
                       require_contraction=True)
service = PreparedValidated(packet)
# service.simulate(initial, inputs), with inputs [20, 1, scenarios]
```

For supplied matrices without numerical discovery, use CLI kind `validated_box`; `examples/validated_box.json` compiles and replays using only the standard library. Replay any compressed external packet with `gzip.open(path, 'rt')`, JSON parsing and `replay_validated`.

## Validation and limits

All 159 tests pass with numerical extras; 142 pass and seventeen optional tests skip in a standard-library-only run. Twenty-five randomized rational models compare the new brackets with exact sharp bounds and execute their retained corners exactly. An independent vertex enumeration checks an arbitrary encoder with two outputs, including a nonzero center and unequal radii. Tests cover precision refinement, source reuse, zero horizon, exact-grid arithmetic, tampered witnesses and statuses, incompatible caches, contraction refusal, nonnormal numerical proposals, nested metric reuse, isolated prepared packets and exact endpoint membership. The CLI example passes the existing compile-to-JSON-to-replay test route.

The next limits are concrete: long-horizon interval widening, spatially uncertain or noncommuting switching, continuous-time intersample accuracy, physically calibrated operating boxes and customer workloads remain outside this new acceptance route. A robustly stable reduced model still needs an independently accepted error contract. The v0.7 all-time engine retains its exact observable-quotient prerequisite, which can be expensive or full-sized on generic dense external models; 0.8 does not remove that limitation. A validated dyadic replacement for its generic full-model tail analysis is a distinct future task. No universal nonlinear reduction, industrial deployment or broad commercial performance lead is claimed.
