# DynaComp: discovering sufficient state for evolving computations

**Complete edition 0.9.0, 10 October 2026 UTC.** This monograph describes an executable target-driven dynamics compiler, from its original rational closure machinery through the current infinite-source extension. It covers exact linear/sparse and polynomial observable closure, finite-machine and weighted-automaton reduction, minimal realizations, parameter loci, finite-horizon uniform and switching bounds, exact stability checks, all-time and correlated-domain bounds, external engineering benchmarks, stable numerical discovery and inexpensive outward-dyadic acceptance. The new spectral chapter accounts for omitted infinite modes before a supplied finite model reaches the reduction engine. Historical chapter versions are retained as such. Claims always concern their stated source formulas, output domains and time horizons; physical calibration and universal nonlinear reduction remain open.

## The question behind the compiler

A large state representation is often designed to describe a system, not to answer one particular question economically. A thermal simulation may represent every cell even when the requested output is an average temperature. A discrete controller may retain implementation details that never affect an observable response. A recurrence may be replayed one step at a time even when its repeated update admits a compact compositional description. The compiler asks which distinctions among states influence the specified future outputs. It seeks an evolving representation that keeps those distinctions and discards the rest.

This question differs from reconstructing the full original state. An output-preserving reduction can remove an unstable hidden mode if the mode never affects the requested output. Conversely, a small hidden variable can be indispensable if it later drives the output. A representation that predicts the current observation perfectly may therefore fail immediately at the next step. The target, admissible updates, initial-state domain and time horizon are part of the reduction problem, rather than incidental metadata attached after compression.

The older IMP and Arithmetic Shadow of Computation ideas supplied the ambition of avoiding unnecessary repeated evolution. Their early projection shortcut did not provide a general solution. Near-unit singular values do not generally identify invariant modes, and an orthogonal projection onto right singular vectors need not equal the limit of repeated application of an operator. The new implementation retains the ambition while replacing that shortcut with explicit closure identities and finite-horizon residual analysis. Its tests include periodic evolution, Jordan growth, oblique fixed-point behavior and stable nonnormal transient amplification.

## Exact linear observable dynamics

Consider a discrete-time rational system with state update x_next = Ax + Bu and output y = Cx. The proposed encoded state is z = Ex. An exact reduced system has update z_next = Gz + Hu and output y = Dz. Three identities establish the relationship: EA = GE, EB = H and DE = C. Substituting the original update gives E(Ax + Bu) = GEx + Hu. The output identity gives Cx = DEx. Induction then establishes matching outputs for every initial state and every finite input sequence. The claim concerns the supplied matrices and target; it does not establish that the supplied model accurately describes a physical installation.

The compiler starts from the row space of C and repeatedly adds its pullback under A. The resulting sequence is the span of the rows of C, CA, CA² and subsequent powers. Dimension increases can happen at most n times for an n-dimensional state space. Once a step adds no new direction, the row space is invariant under right multiplication by A. A rational row basis becomes E, and exact coordinate solving supplies G and D. H is simply EB. All rank decisions and coefficient calculations use Fraction arithmetic, without a floating-point rank threshold.

This construction is minimal among linear encoded states that preserve the requested outputs for all initial states and all input sequences. Any such encoder must contain the rows of C in its row space because the current output factors through it. Transition closure then forces that row space to contain CA and every subsequent pullback. The compiler therefore constructs a space contained in every admissible linear observable representation. Minimality does not imply that nonlinear encodings, restricted initial-state sets or approximate models cannot be smaller.

Its kernel also has a behavioral interpretation. Two initial states with difference v produce identical outputs under the same future inputs exactly when CA^t v vanishes for every t. Input contributions cancel in the difference. The Cayley–Hamilton identity reduces this all-future test to the first n powers. The compiler's encoder contains precisely the row space generated by those powers, so Ev = 0 expresses the same indistinguishability. The repository provides an earliest-distinguishing-horizon calculation over those first n observations.

## Failure witnesses and representation repair

Suppose an initial encoder E does not close under A. There is a vector v in its kernel with EAv nonzero. The states zero and v have the same encoded state now but different encoded states after one update with the same input. No deterministic update function of the current encoding alone can match both. This is a useful failure explanation: the representation discarded a distinction that its next update needs.

The linear compiler records such a witness before each dimension-increasing refinement. It then enlarges the row space with the missing pullbacks. In the hidden-driver example, the target is the first coordinate, its update depends on a second coordinate, and a third coordinate evolves independently with multiplier two. The compiler keeps the target and its driver, removes the irrelevant third coordinate, and records why the target alone was insufficient. The repair is exact for a known rational linear model. It is not yet an adaptive algorithm for arbitrary approximate nonlinear lifts learned from trajectories.

Independent replay of a linear packet reads the original and reduced matrices and checks the three identities directly. It also checks that the encoder has independent rows. Discovery can therefore change without changing the basic acceptance contract. Replay uses Python's rational arithmetic and ordinary program execution. No Lean theorem, proof-kernel acceptance or verified implementation of Fraction arithmetic is asserted in this release.

## Sparse interactions and discovered aggregates

Dense storage is unnecessary when the update operator is sparse. The sparse backend stores rational nonzero triples and computes observable pullbacks directly from them. Its exact acceptance identities have the same meaning as the dense backend. A user-selected observable rank limit prevents an unexpectedly large closure from silently becoming an expensive dense algebra problem. Exact rational arithmetic remains best suited to moderate observable ranks; sparse input alone does not make every symbolic calculation scalable.

A complementary discovery routine searches for closed group sums. Start with a partition of coordinates that respects the requested target groups. For each coordinate, compute how its value contributes after the update to every current group sum. Coordinates can remain in the same class only when these contribution vectors agree. Repeated refinement reaches the coarsest closed aggregate partition refining the initial partition. This is a structural aggregation statement, not a general graph-isomorphism or symmetry-group classification algorithm.

The regular periodic ring demonstrates the simplest closure. Its diagonal weight is 49/100 and its two neighboring weights are each 1/4. Every column sums to 99/100, so the total state sum evolves by that scalar multiplier. A 64-coordinate example therefore reduces to one observable state. A second example inserts a local diagonal defect into a 16-coordinate symmetric ring. Aggregate refinement then discovers nine groups, retaining the reflection pairs around the defect. The representation changes in response to the interaction structure rather than assuming that the old uniform aggregate remains valid.

## Nonlinear updates with exact feature closure

For an autonomous polynomial update F, an observable polynomial p evolves by pullback p composed with F. A finite collection of polynomial features closes when every pulled-back feature lies in their linear span over the rationals. The polynomial compiler starts from the requested target polynomials, expands their pullbacks, and performs exact row-space reduction in a monomial coefficient representation. If the span stabilizes, it emits a linear feature update and a linear target decoder.

For x_next = x + y² and y_next = y, the target x generates the feature y². The two features x and y² then evolve linearly. The sign of y is unnecessary for this target, even though reconstructing the full original state would need it. This example demonstrates a change of sufficient coordinates, rather than a reduction from two coordinates to one. The benefit of a lifted representation can include simpler update operations or reuse across many target queries, and must be measured rather than inferred from a feature count alone.

The compiler has explicit degree, feature and round limits. An update such as x_next = x² repeatedly doubles the degree of the target x and reaches the configured expansion limit. The returned status states that the search exceeded resources. It does not prove that no alternative finite representation exists. Polynomial identity replay checks the emitted pullback and output relations coefficient by coefficient; finite trajectory agreement is not used as a substitute for those identities. Controlled polynomial systems and domain-specific relations among variables are not implemented in this version.

## Finite machines and future equivalence

The finite backend accepts a deterministic partial transition table and an observation for every state. An unavailable action is itself observable: it is distinct from taking a defined transition whose destination happens to have an ordinary observation label. The initial partition groups equal observations. Refinement splits states according to successor classes and transition definedness until the partition stabilizes. The quotient then preserves the observation and every available action.

Preservation establishes that states in one quotient class have the same future traces. To establish that distinct retained classes genuinely differ, the packet contains a distinguishing action word for every pair of classes. Breadth-first search supplies a shortest word. Independent replay checks the quotient homomorphism and executes each pair-specific word. This handles larger machines where one experiment cannot distinguish every pair, rather than selecting one witness from a small fixture and pretending it is universal.

The regression suite exhaustively processes the 324 two-state machines with binary observations and two partial actions. Each of the four transition slots can be unavailable or point to either state, yielding 81 transition tables, and there are four observation assignments. The suite also includes a larger machine with duplicated equivalent states and different distinguishing experiments. These checks concern explicitly supplied finite models. They do not establish the behavior of a black-box implementation outside its observed or modeled transition space.

## Approximate reduction and finite-horizon error

An approximate model evolves z_next = Gz + Hu and reconstructs a candidate full state Vz. Define the error e = x − Vz, transition residual R = AV − VG and input residual S = B − VH. Direct substitution gives e_next = Ae + Rz + Su. Repeated substitution gives e_t = A^t e_0 plus the sum of propagated residual terms A^(t−1−j)(Rz_j + Su_j). Multiplication by C turns this into the output-error identity.

The implementation bounds each propagated residual term using the induced infinity norm and the exact rational magnitude of its reduced state or input. The initial term CA^t e_0 is evaluated exactly before taking its norm. These bounds include hidden initial-state error, which a small transition residual alone cannot eliminate. They apply to a supplied finite input trace and initial state. They do not automatically cover every initial state in a population, unknown future input ranges or parameter uncertainty.

The method does not require asymptotic stability. The nonnormal example A = [[1/2, 8], [0, 1/2]] has decaying eigenvalues, yet an initially hidden second coordinate immediately contributes eight units to the observed first coordinate. The output envelope correctly records that first-step error. A requested tolerance of one therefore gives a usable prefix ending at the initial time. Later recovery below tolerance does not retroactively make the intervening failed interval usable.

Positive weights also define a weighted infinity norm. The exact factor is the maximum weighted absolute row sum. If the factor is below one, the matrix is contractive in that norm. A simple resolvent-like search proposes weights and accepts only their rationally calculated factor. This can establish a margin for the nonnormal example even though its unweighted row norm exceeds one. Failure of diagonal metric search is not a proof of instability: a system may need a more general metric, a different norm or a different analysis.

## Numerical discovery and compositional fast-forwarding

Snapshot-based POD is an optional proposal mechanism. It computes numerical singular vectors, rationalizes the proposed lift, checks independence and constructs an exact rational left inverse. The resulting reduced matrices are labeled as a candidate requiring error analysis. Singular values on the training snapshots do not prove global closure, and rationalization can alter the candidate geometry. This distinction keeps numerical discovery useful without transferring an unsupported guarantee from a data fit to every future state.

Repeated affine updates offer a separate acceleration. A fixed input word induces an affine map on the reduced state. Composition of two maps is exact, and binary exponentiation applies that block a large number of times with logarithmically many matrix compositions. The trillion-step integrator example returns the exact target 500000000002500000000002 from initial state (2,3) and unit driver input. Large rational numerators and denominators can still make arithmetic expensive; a logarithmic composition count is not a constant-time or bounded-bit-complexity claim.

## Evidence, performance and the next research step

The current implementation has 131 test methods, including the network and thermal applications and optional numerical backend regressions. Earlier release evidence is retained below, followed by the merged extension chapters. Beyond the exhaustive finite corpus, these include 60 random rational forced-system comparisons, 30 random output-envelope comparisons, 24 random forced-memory models, uniform-box corner comparisons, adaptive rank completion, switching-family trajectories, mutation rejection for emitted identities, nonlinear feature trajectories, sparse aggregation, contraction margins and affine-block replay. Deterministic examples regenerate retained receipts. Tests support the implementation of the stated mechanisms, while the mathematical derivations explain their scope and assumptions.

The optional periodic-stencil benchmark uses 2,048 coordinates, 256 initial-state queries and 200 steps. The retained local receipt reports about 0.478 seconds for full batched stencil evolution, about 0.000422 seconds for reduced evolution including initial encoding, and about 0.121 seconds for sparse construction and identity replay. Including construction, the first batch is about 3.94 times faster. The maximum floating output difference across every sampled time is about 5.68 × 10^−14. The full and reduced paths are warmed and timed in alternating order, with three measurements and median timings retained.

This performance result is intentionally narrow. It preserves an exactly closed aggregate of a particularly structured model, and does not reconstruct the field. It does not establish an advantage over established model-reduction packages or over a hand-written aggregate recurrence. Its significance is that the compiler discovers and checks the reduction from the supplied sparse update, and the measurement includes its construction cost. Future evaluations need harder target sets, defects, forcing, parameter variation and matched-tolerance baselines.

The continuation implements a first cost-aware enrichment route and distinguishes structural family closure from accidental specialization. The strongest next research direction is to extend that route beyond coordinate selection, determine when constructing a repaired model pays off on repeated queries, and handle nonlinear memory. A physical thermal or lattice application can make those tradeoffs measurable. The eventual ambition is a compiler that discovers sufficient evolving state and knows the domain and horizon of its replacement; the present repository supplies a precise, executable foundation for that effort.

## Uniform guarantees for populations and future inputs

The first release bounded the error of a specified initial state and input trace. A useful reduction decision often needs to cover a population of initial states and future controls that have not yet been chosen. Version 0.2 supplies this extension for rational linear models. The initial state lies in a componentwise box centered at a supplied vector, and each future input has a supplied nonnegative componentwise magnitude bound. The reduced initial state is Ex, and the lift and encoder satisfy EV = I. The guarantee covers every member of that box and every input word satisfying the declared bounds.

The initial reconstruction error is (I−VE)x. Its central component is evaluated exactly, while its uncertain component is bounded by the absolute matrix applied to the initial radii. Reduced-state magnitude bounds propagate by |G| and |H|. The output-error identity from the earlier section then uses absolute propagated residual matrices applied to those componentwise state and input envelopes. This is a conservative interval-style inequality. It does not assume independence among original coordinates, and it can overestimate error by ignoring correlations. Conservatism preserves coverage but may require more retained coordinates than a sharper analysis.

The implementation stores CA^k directly rather than full A^k matrices for this route. That reduces one source of work when the output dimension is small. The calculation still uses dense rational algebra and a finite-horizon convolution, so the current method is intended for modest model dimensions and horizons. It is not a scalable sparse interval solver. The regression suite checks the bounds against every initial and input corner in a small multi-output example, as well as a noncoordinate rational lift with a nonzero center.

## Adaptive enrichment and its economic objective

For a chosen coordinate subset, the encoder selects those coordinates and the lift inserts them into the full space with zeros elsewhere. The candidate dynamics are the corresponding projected operators. DynaComp calculates their uniform output envelope over the requested horizon. If the largest bound exceeds tolerance, it tries each missing coordinate, measures the reduction in that largest bound and divides by the supplied coordinate cost. The greedy route selects the best score, reports its trials and continues until tolerance is met or the rank budget is exhausted.

The driver demonstration begins with one retained coordinate and a worst bound of 13/4. Adding the hidden driver reduces that envelope to zero. A third coordinate has unstable multiplier seven, but cannot affect the output and is not retained. The eight-state attenuated-transport demonstration shows a less trivial tolerance tradeoff: ranks two, three, four and eight satisfy tolerances 1/10, 1/100, 1/1000 and zero, respectively, over eight steps and the declared input radii. These are bounds covering a population, rather than errors observed on a few favorable trajectories.

The greedy score does not establish optimality. A retained counterexample makes this visible: greedy enrichment incurs additive cost seven, whereas another accepted coordinate subset costs six. The best-first route enumerates coordinate subsets in increasing additive cost. Because costs are strictly positive, the first subset passing the envelope criterion has least cost among the permitted subsets, with the declared required coordinates and rank limit. No monotonicity of error under enrichment is assumed. An evaluation limit yields an incomplete-search status, and a completed unsuccessful search says only that no subset passed this conservative criterion. Neither status establishes impossibility for the underlying dynamics or for other kinds of representation.

## Structural families rather than accidental parameter values

An exact quotient at one coefficient choice can fail as soon as a coupling changes. DynaComp now closes the target row space under every coefficient operator in a supplied family. The resulting encoder obeys EA_i = G_iE and EB_i = H_i for each coefficient. Linearity then gives the corresponding identities for arbitrary coefficient combinations. The same encoder remains valid even when those coefficients vary from one discrete time step to the next. This establishes exact output preservation, not stability of every coefficient combination.

A second routine takes a fixed encoder and an affine operator family A0 plus parameter-weighted coefficient matrices. It constructs an exact right inverse and calculates each coefficient's closure defect. Every scalar entry of the resulting affine defect must vanish. The returned parameter equations therefore describe where that fixed representation closes, without inferring the condition from sampled parameter values. In the retained example, the one-state encoding closes only when the coupling parameter is zero, while a two-state encoding closes across the family. General determinantal rank stratification and nonlinear parameter dependence remain outside this version.

## Output memory as a derived mathematical object

An output can have a finite exact recurrence even when its current value alone does not determine its next value. For a nonzero scalar linear target, DynaComp forms the rows C, CA and successive pullbacks until the next row is a linear combination of its predecessors. Their first dependence determines the minimal output recurrence order over all initial states. Expanding each output under future forcing also determines input-memory kernels. These coefficients are obtained from the supplied operators and are checked by rational identities, rather than estimated from a finite time series.

For the forced integrator, the result is y[t+2] = 2y[t+1] − y[t] + u[t]. The second input kernel is zero because the latest control changes the hidden driver before it affects the target. The recurrence simulator requires its initial output window to be consistent with the chosen initial input prefix. Independent replay checks the state-elimination identity, the input kernels and independence of the initial observable rows. The random regression corpus compares these recurrences with full evolution for 24 models with two input channels. This derives linear memory exactly; discovering a sufficient finite memory representation for arbitrary nonlinear projections remains an open development problem.

## An operational application: repeated inventory-network scenarios

Version 0.3 turns the aggregate machinery into a focused scenario application. Its state is stock distributed over nodes of a fixed routing network. Its inputs are replenishment amounts sent through named channels with explicit distributions over nodes. Its outputs are stock totals and dispatch quantities for each reporting region. A scenario supplies an arbitrary initial stock vector and a time-ordered replenishment sequence. The application returns every-time regional outputs, cumulative dispatch, minimum margins above supplied stock floors and violation counts. These are operational quantities a user can compare across scenarios, rather than only matrices describing a reduction.

The model has a nonnegative, column-substochastic transition operator A. Every source node's outgoing fractions sum to at most one. The missing fraction is interpreted as dispatch from that source region. This interpretation must match the supplied physical model: the same mathematical exit could otherwise mean leakage or spoilage. Each input channel has a nonnegative node distribution summing to one, so replenishment quantities retain their mass meaning. At time t stock is observed, routed and dispatched before input u[t] arrives. A horizon of T inputs therefore produces T+1 stock observations and T dispatch observations.

Let the missing transition fraction at node j be ell[j] = 1 minus the sum of column j. The full system obeys x[t+1] = A x[t] + B u[t]. Total stock obeys the identity sum(x[t+1]) = sum(x[t]) minus sum(ell[j]x[t,j]) plus sum(u[t]). Telescoping gives exact conservation: initial stock plus all replenishment equals final stock plus cumulative dispatch. The report calculates this residual. Rational execution makes it zero; ordinary floating-point execution can introduce rounding error. The numerical report does not turn this algebraic identity into a certified floating-point enclosure.

The initial aggregate labels record each node's reporting region and exact dispatch fraction. Refinement compares each source node's total contribution to each current destination class. Two nodes remain together only if these contribution vectors agree and their earlier labels agree. Repeating the process yields a closed aggregate partition. Heterogeneous internal routing is allowed: nodes can route differently inside a class while contributing the same total to each class. A defect that changes aggregate contributions splits the relevant class and can propagate further refinement. A group budget stops overly expensive discovery and returns a resource-limit exception, without accepting an approximate quotient or claiming that all other representation families have been exhausted.

The compact plan stores one class index per node, the source model and rational G, H and D. It avoids serializing the dense binary encoder E. Independent replay reconstructs class-column sums and checks region/dispatch consistency, then recomputes reduced input and output coefficients. These conditions establish EA = GE, EB = H and DE = C. Induction proves identical regional stock and dispatch at every future time for every initial state and input word. Positivity and replenishment normalization supply the conservation interpretation. Discovery produces the coarsest refinement of the declared labels, but replay establishes output preservation rather than independently proving minimality of an asserted refinement history.

The application has a dependency-free scalar backend for rational or binary64 execution and an optional NumPy/SciPy backend for scenario batches. The numerical backend groups scenarios by horizon, preserves their original order and uses sparse original evolution or compiled evolution as requested. The benchmark uses this same production executor. Stock-floor checks are diagnostic and do not influence routing. A floor comparison near equality can differ through rounding, so exact execution is the appropriate route when the decision margin is small. Capacities, congestion, demand caps, backlog and empirically learned transition fractions are not implemented here.

## Performance with construction costs and an explicit losing case

The retained network benchmark is synthetic. It contains four reporting regions with heterogeneous internal routing but common regional transfer and dispatch fractions. The structured cases have 128, 2,048 and 8,192 original states, all compiling to four states. Each evaluates 128 scenarios through 96 steps. Full evolution uses SciPy compressed sparse multiplication, not Python iteration, and records all regional stock and dispatch outputs. An independent hand-built four-compartment model supplies a second baseline. Five warmed measurements rotate execution order, and the receipt retains individual timings and medians.

On the retained machine, the 8,192-state case takes about 0.296 seconds for a full sparse batch and 0.00143 seconds for compiled execution, including initial encoding and every-time outputs. Rational construction and identity replay take about 0.714 seconds. Including numerical setup and construction, the first compiled batch is slower: its speed ratio is about 0.784. The estimated amortization point is about 195 queries in similarly sized batches. The hand-built aggregate takes about 0.000885 seconds, faster than the automatic executor. The largest full/compiled output difference across sampled times is about 2.55 times 10 to the minus 10. Exact output equivalence comes from rational identities; this observed numerical difference measures ordinary binary64 execution.

A deliberately heterogeneous 32-state case assigns unique per-cell dispatch rates. Its requested outputs prevent aggregation, so it retains all 32 states. Compiled execution is slightly slower and incurs additional construction cost. This case matters: the application identifies a mathematically valid full-size representation rather than manufacturing a compression benefit. Neither successful structured cases nor this losing case establishes a result on customer data. Model generation, JSON I/O, scenario ingestion and dependency installation are excluded from the kernel benchmark, and a customer evaluation must include the costs relevant to its actual workflow.

The commercial hypothesis is repeated analysis of existing stock-flow models whose requested outputs admit substantial closure. A paid pilot should import a representative customer model, preserve agreed outputs, reproduce its baseline and show an improvement over an agreed workload after accounting for all integration and construction costs. Where a customer already has a correct aggregate, the value must come from automation, coverage across models or maintenance rather than an invented numerical advantage. The repository now provides a concrete application and a reproducible evaluation route. Customer fit, measured industrial performance and recurring revenue remain to be established.


# DynaComp 0.4: realization, automata, parameter loci and direction enrichment

**Companion chapter to the [monograph](MONOGRAPH.md), 10 October 2026 UTC.** Version 0.1 compiled the observable quotient of a supplied model. Versions 0.2 and 0.3 added box-uniform adaptive coordinates, coefficient families, input memory and an inventory-network application. This release adds five more mechanisms. Each one keeps the same contract: a compiler discovers a reduction, and an independent replay accepts it by checking exact identities. Each section below gives its derivation, its acceptance identities and its limits.

## Minimal realization and the Hankel lower bound

The observable quotient keeps every distinction that *some initial state* could make. If the system starts at rest, a state that no input sequence can reach never influences the output. Such a state is redundant even when the target observes it. `minimal_realization` therefore first takes the reachable column space. A breadth-first Krylov search over the columns of B, AB, A²B, … gives a basis P with **AP = PF** and **B = PK**. The compiler then applies observable closure to (F, K, CP), which gives E, G, H and D with **EF = GE, EK = H and DE = CP**. Substitution gives CAᵏB = CPFᵏK = DEFᵏK = DGᵏEK = DGᵏH for every k. The Markov parameters agree, and so does the zero-initial-state output for every input sequence.

Minimality is certified independently of how the quotient was found. Any realization of order s factors the block Hankel matrix [CA^(i+j)B] through an s-dimensional space. Its rank is therefore at most s. The packet names s output-power pairs (i, k) and s input-power pairs (j, l). Replay recomputes the s × s block of entries C_i A^(k+l) B_j from the *original* matrices and checks that the block is nonsingular. A nonsingular block shows that no realization can be smaller than s, and the identities show that the emitted one has order exactly s. The compiler chooses these experiments by running the same labeled Krylov search on the minimal system. Their observability and reachability blocks are both nonsingular, so their product is nonsingular too.

Two examples show the gain. In `examples/excitation.json`, observable closure keeps three states. Only two of them are reachable from the actuator, so the input-output map has order two. In a separate test, a diagonal model observes both of its modes, but its actuator excites only one of them, so the input-output map has order one.

## Switched systems and weighted automata

A switched system applies one of several updates A_s at each step. To be exact for every mode sequence, the encoder's row space must contain C and be invariant under every A_s at once. `compile_switched` builds this space by breadth-first search over words. Each retained row is c_i A_w for a recorded word w, so the encoder is also a set of experiments: "observe output i after applying modes w". Replay checks E A_s = G_s E and E B_s = H_s for each mode, and DE = C. This is the same span that `compile_family` computes for a coefficient family, and a randomized test checks that the two agree. The switched form adds per-mode inputs, the experiment labels, and the automaton and equivalence machinery that follow. Sequential closure under one mode at a time is not enough. Under the identity mode alone, the target e₁ needs one state. Adding a mode that swaps the two coordinates makes a second state necessary.

A weighted automaton (α, {A_s}, η) assigns each word w = s₁…s_k the weight αA_{s₁}⋯A_{s_k}η. Probabilistic automata, hidden Markov models (with A_s[i][j] = T[i][j]·O[j][s]) and rational power series are all special cases. `minimize_automaton` takes the forward span of the vectors αA_w and then the backward span of the vectors A_wη inside that quotient. This is Schützenberger's minimization over Q. The certificate has two parts. The forward and backward identities show that the small automaton assigns every word the same weight as the original. A nonsingular Hankel block of weights f(u_i v_j), over recorded prefix and suffix words, proves that the result is minimal among all weighted automata over Q.

`check_equivalence` decides whether two automata define the same weighted language. It explores the joint span of the vectors (α₁A₁_w, α₂A₂_w) in shortlex order, and evaluates the weight difference on every word it visits. The first word with a nonzero difference is a *shortest* distinguishing word. The argument is as follows. Let L be the length of a shortest distinguishing word. Every word's vector lies in the span of retained words no longer than itself. So some word of length at most L, obtained by extending a retained word by one letter, also distinguishes, and that word is visited before any longer word. If no visited word distinguishes, the explored span is invariant and annihilates the difference vector (η₁, −η₂). The packet records that span as an equivalence certificate, and replay checks its invariance under every letter. The randomized tests compare both outcomes against brute-force enumeration of words.

In `examples/hmm_equivalence.json`, a three-state weather HMM splits its Sunny state into two identical copies. It is certified equivalent to the two-state model over all observation sequences of every length. Changing one copy's emission probabilities by 1/50 breaks equivalence, and the packet names a one-symbol word on which the two models disagree.

## Structural versus accidental closure

The development notes identified parameter dependence as the next substantial step. A reduction may close at one coefficient choice and fail nearby. `closure_conditions` (0.2) derives affine equations for where one *fixed* encoder closes. `parametric_analysis` accepts A(t), B(t) and C(t) whose entries are polynomials in one real parameter t. It analyzes three matrices: the observability stack O(t), the reachability stack R(t) and the Hankel product O(t)R(t), whose rank is the minimal order.

For a real polynomial matrix M(t), let e_k(t) be the k-th elementary symmetric function of the eigenvalues of M(t)ᵀM(t). By the Cauchy–Binet formula, e_k is the sum of the squares of all k × k minors of M. For real t, rank M(t) ≥ k therefore holds exactly when e_k(t) ≠ 0. The generic rank r is the largest k for which e_k is not the zero polynomial. The exceptional real parameters are exactly the real roots of e_r. The implementation computes every e_k by Faddeev–LeVerrier recursion over Q[t], with no evaluation at sample points. It then reduces e_r to its square-free primitive part, and counts and isolates the real roots exactly with a Sturm sequence. A rational root N/L of an integer polynomial with leading coefficient L makes Lx an integer. Refining each isolating interval below width 1/L therefore leaves at most one integer candidate, so exact rational exceptions are found without factoring any coefficient. Irrational exceptions are reported as rational isolating intervals.

In `examples/symmetric_rod.json`, a three-cell rod has its heater and sensor in the middle cell. When the far end's self-coupling t equals the near end's value 1/2, the rod is mirror-symmetric. The antisymmetric mode then disappears from the input-output map, and the minimal order drops from three to two. The analysis returns t = 1/2 as the only real exception and labels a reduction compiled there as **accidental**. This is how symmetry candidates should enter: they appear as checkable rank drops of the actual operator, not as assumed structure.

Replay re-derives each profile from the stored model. It then checks, by direct rational elimination at specific parameter values, that the rank at the recorded generic witness equals the generic rank and that each rational exception lowers it. The analysis currently supports one parameter. Several parameters need multivariate elimination and remain open.

## Certified greedy enrichment

The monograph named target-weighted adaptive enrichment as the strongest next research direction. `adaptive_coordinates` (0.2) adds original coordinates with guarantees that hold uniformly over a box. `greedy_enrichment` goes beyond coordinate subsets: it adds general rational directions, with guarantees for each supplied scenario. It starts from the target directions (the rows of C taken as columns of V). It uses the exact orthogonal projector W = (VᵀV)⁻¹Vᵀ, with G = WAV and H = WB. For every supplied scenario (an initial state and an input trace), it computes the rational finite-horizon envelope from `bounds.finite_horizon_bound`. While the worst envelope exceeds the tolerance, it adds the full-state error snapshot x_t − Vz_t whose future is most visible to the target, scored by Σ_k |CAᵏe| over a lookahead. If every snapshot error vanishes while an envelope is still positive, it adds a column of the transition residual AV − VG instead.

Soundness does not depend on how a direction is chosen. Every returned envelope is recomputed exactly for whatever rational V results. The algorithm can therefore round each new direction to small denominators to keep rational growth under control without weakening any certificate. Once V spans the whole space, every residual and the initial error vanish. The loop therefore always ends, with status "certified" or, if a dimension budget is reached first, "budget_exhausted". Replay checks that WV = I, that G and H are the stated projections, and that every stored envelope matches a fresh computation.

The worst certified bound is not monotone in the dimension. A Galerkin projection can make the transient worse before it makes it better, and the receipts keep that history rather than hiding it. The certificates cover the supplied scenarios only. Covering a population of initial states, a set of input ranges or parameter uncertainty needs a different bound.

## Controlled polynomial closure

Version 0.1 handled autonomous polynomial updates only. `compile_controlled_polynomial` accepts x_next = F(x, u), where u is an arbitrary input at each step. It pulls each feature back through F and splits the result by input monomial: p∘F = Σ_m u^m c_m(x). Because the inputs are arbitrary, every coefficient polynomial c_m must lie in the span of the features. The closure loop adds them until the span stabilizes. The result is an exact lifted system z_next = Σ_m u^m G_m z, which is linear in the features and polynomial in the inputs. Replay recomputes every pullback and compares it coefficient by coefficient.

For x_next = x + u·y² and y_next = y/2 + v, the target x closes with the four features x, y², y and 1. The input v enters y² quadratically, so the lift carries blocks for 1, v, v² and u. The same degree, feature and round limits apply as in the autonomous case. Reaching a limit does not show that no finite closure exists.

## Engineering changes

The sparse compiler now uses an incremental echelon form. It tests each new pullback once, instead of re-reducing the whole span every round, and it reads coordinates in the final reduced basis directly from its pivot columns. On a 16 × 16 heat stencil this cut compile time from about 38 s to about 1.3 s on the development machine, with identical packets. Every previously retained receipt, including the network application's, regenerates byte for byte. Computing Hankel blocks incrementally made realization replay about three times faster.

The thermal plate application in `examples/thermal_plate.py` compiles a 12 × 12 heated plate whose target is the summed temperature of the four corner cells. The observable order is 17 of 144 states, the closed aggregate partition has 21 classes, and the Hankel-certified minimal order is also 17. Reduced and full evolution agree *exactly*, in rational arithmetic, for 40 steps of pulsed heating from an arbitrary initial field. The reduction preserves the corner sum only; it does not reconstruct the field.

## What is still not claimed

No multivariate parameter analysis, no nonlinear memory closure, no continuous-time generators and no general Lyapunov metrics are implemented. Every claim of minimality, equivalence or exceptionality concerns the supplied model over Q, or for parameters over the real numbers. None of these claims extends to a physical installation the model is meant to describe. Replay is ordinary Python execution; no proof-assistant acceptance is claimed.


# DynaComp 0.5: sharp uniform errors, general directions and sparse uncertain families

Version 0.4 supplied general-direction enrichment with guarantees for retained scenarios. Earlier box-uniform acceptance was attached to coordinate selection. Version 0.5 combines these capabilities and strengthens the acceptance calculation: for a fixed rational linear model, it calculates sharp error suprema over independently chosen componentwise initial-state and input boxes. Sparse target pullbacks keep full matrix powers out of this route. A separate sparse family compiler extends the switched/common-encoding machinery to supplied sparse operators. The thermal workbench demonstrates both mechanisms without conflating their domains.

## Direct output-error kernels

Write the full system as x[t+1] = A x[t] + B u[t], y[t] = C x[t]. A candidate has z[0] = E x[0], z[t+1] = G z[t] + H u[t], and output D z[t], where D = C V and EV = I. The output error at t is exactly K[t] x[0] plus the sum over j < t of J[t−1−j] u[j], with K[t] = C A^t − D G^t E and J[k] = C A^k B − D G^k H. This follows by expanding both recurrences, and it does not require stability, approximate closure or any probabilistic assumption.

Let the initial state be c plus componentwise deviations bounded by radius r. Let input j have center a[j] and nonnegative radius b[j], with deviations chosen independently across coordinates and times. For output i at time t, define mu as row i of K[t] applied to c plus the sum of row i of J[t−1−j] applied to a[j]. Its exact maximal absolute error is |mu| plus the sum of |K[t,i,k]| r[k] plus the corresponding absolute input-coefficient sums. All quantities are rational. The result is sharp because each free deviation can take the sign of its coefficient, oriented to reinforce mu. If mu is zero, either common orientation works.

The packet retains a witness at the worst output/time: one initial state and one full input word within the declared boxes. Its signs attain the stored maximal absolute error. Replay reconstructs the kernels and domain, recomputes all componentwise bounds and the witness, then independently simulates that witness through the full and reduced systems. It checks the signed error and attainment. An above-tolerance bound therefore identifies an actual failing scenario for this fixed candidate, rather than only a conservative envelope. It does not establish that every candidate of the same rank fails.

Sharpness applies to each component/time range and to the maximum over the finite horizon. Different output/time extrema can require different input words, so the table is not a collection of independently attainable values along one trajectory. Correlated initial-state constraints, coupled controls, uncertain coefficients and nonlinear dynamics require different calculations. Binary64 arithmetic and QR rounding are also outside this exact rational error guarantee.

## Why the direct form matters for general directions

Residual-propagation envelopes can discard cancellations before the target sees them. Consider A = one-half times the identity, target C = (1,1), encoder E = (1,1) and lift V = (1/2,1/2) transposed. This one-state representation preserves the target exactly. Yet first taking absolute values of the hidden reconstruction map and then applying the target can create a positive initial-error envelope. The direct kernel is C(I−VE) = 0, so the new bound correctly returns zero at every time. The retained regression exposes this difference, rather than assuming all conservative envelopes behave well under noncoordinate projections.

The implementation caches C A^k using sparse row pullbacks. It evolves D G^k in the small output/reduced space and forms the two error kernels. It never constructs A^k. Input convolution still costs quadratically in the horizon in the current implementation, and rational coefficients can grow. These changes remove a major dense-state cost but do not make arbitrary large horizons inexpensive.

The existing scenario-specific residual bound also now evolves C A^k directly and caches propagated residual norms. Its mathematical envelope and emitted receipts are unchanged. This improves the machinery used by the scenario enrichment module without changing its domain or silently replacing its conservative formula with a different criterion.

## General-direction enrichment with uniform acceptance

`uniform_enrichment` starts from the target row space and optional supplied seed directions. With full-row-rank E, it uses the rational orthogonal right inverse V = Eᵀ(EEᵀ)⁻¹, sets G = EAV, H = EB and D = CV, and evaluates the sharp box certificate. Consequently its state directions need not be original coordinate axes. Initial states and inputs are not training samples: every admissible choice in the stated boxes is covered.

The worst-kernel strategy adds an independent missing output sensitivity from the worst time/output. If that row is already dependent, it falls back to the earliest missing target pullback. The prefix strategy always adds the earliest missing pullback. Independence is checked incrementally. These are heuristics, and the sharp worst-case error need not decrease monotonically when the projection changes. The history records each evaluated rank and bound. A dimension budget reports `rank_budget_exhausted`; a successful candidate reports `tolerance_met`.

For a fixed finite horizon, containing every target pullback through that horizon suffices for zero error, even if the row space is not invariant for all future times. This supplies finite termination when the budget permits the needed span. A retained three-state shift example closes through one step with two states, while the all-future target requires three. The certificate explicitly records its horizon, so the shorter guarantee cannot be silently extended to longer scenarios.

Replay checks the projection coefficients, recomputes the final certificate and verifies status against tolerance and budget. It does not rerun enrichment or accept an optimality claim from the discovery history. Existing least-cost coordinate search remains a separate optimizer within its own representation family and conservative criterion. General-direction minimum-rank or minimum-runtime optimization remains open.

## Common sparse families

`compile_sparse_family` accepts named sparse operator/input pairs with one shared target. A labeled breadth-first Krylov search closes the target row space under every operator, using sparse pullbacks. One final rational reduction gives a common E and D and separate G_s and H_s. Replay checks E A_s = G_s E, E B_s = H_s and DE = C without rerunning discovery.

Linearity implies the same identities for A[t] = sum alpha_s[t] A_s and B[t] = sum alpha_s[t] B_s. They therefore hold under arbitrary time-varying coefficient combinations. When two endpoints have the same input matrix and their coefficients sum to one, their affine interpolation has that same input matrix. This is the interpretation used for thermal conductance. Output preservation does not itself imply positive updates, a stable physical timestep or an accurate physical model. A supplied rank limit constrains discovery, and does not establish that smaller approximate representations are impossible.

Randomized cross-checks compare the sparse common row space with the existing dense switched compiler. A separate exact test uses varying interpolation coefficients, including a coefficient outside the physical interval, to check the algebraic claim. The thermal application restricts coefficients to its physically nonnegative timestep interval.

## Heated-plate workbench

The workbench models excess temperature above a fixed ambient on an even square grid. Neighbor conductance and ambient loss are nonnegative, and 4k plus loss must not exceed one. Four central cells share one heater input, normalized so one unit of input contributes one unit of total excess temperature. The target is the mean of the four corner temperatures. The finite stencil is synthetic; no claim about continuum discretization error or a measured installation follows from its exact algebra.

In the default 12-by-12 case, conductance is 1/8, loss is 1/100, initial temperatures range independently from 18 to 22, and each heater input ranges independently from zero to two for twenty steps. Prefix enrichment produces a nine-state model with a sharp uniform worst output error of approximately 0.0828409143, below the requested 0.1 tolerance. The full model has 144 states. The report also supplies sharp full-model temperature intervals, which describe output uncertainty itself rather than reduction error, and classifies those intervals relative to a supplied threshold. It runs nominal, heater-outage and upper-box-corner scenarios in exact arithmetic.

A separate common quotient has seventeen states and preserves the target exactly when conductance changes arbitrarily at each step within [1/16,1/8]. Its endpoint identities establish this for every such conductance word, every initial field and every heater sequence. The application also simulates a changing interior-conductance word in rational arithmetic. The nine-state approximation and its output intervals fix conductance at 1/8; they do not inherit the larger exact family's parameter guarantee.

The optional numerical batch executor uses a QR basis for the certified subspace to avoid evolving in an ill-scaled rational echelon coordinate system. The public convenience API checks each supplied binary64 value against the rational box endpoints using its exact binary rational value. QR and evolution still introduce uncertified rounding. The benchmark compares this same prepared executor with full SciPy sparse evolution, includes initial encoding and every-time outputs, and separately records construction, numerical preparation and amortization.

On the retained machine, 256 scenarios over twenty steps take about 0.00203 seconds full and 0.000211 seconds reduced, approximately 9.6 times faster in execution. Construction and rational replay take about 0.772 seconds, so the first batch loses substantially and estimated amortization requires about 426 similarly sized batches. The largest observed output difference is about 0.0753, below the rational uniform model-reduction bound. Generation, JSON I/O, ingestion, installation and per-batch domain checks are excluded from the prepared-executor timing. This is a small synthetic kernel comparison, not evidence of industrial profit or a universal performance advantage.

## Validation and remaining work

The complete suite has 117 test methods. New tests compare sharp componentwise bounds with every box corner on sixteen random two-state models, verify output ranges by enumeration, check sparse/dense uniform completion on twelve random models, compare sparse families against twelve random dense switched families, exercise attaining witnesses, reject packet mutations, test finite-horizon closure and budget exhaustion, verify thermal conservation and changing conductance, and check the optional numerical backend and rounded-endpoint domain rejection. Original retained receipts regenerate unchanged except for additive summary entries.

General-direction uniform acceptance is now implemented for fixed rational linear models. Nonlinear memory, multivariate rank stratification, correlated-domain sharp bounds, approximate reductions uniform over uncertain parameters, continuous-time generators, general Lyapunov metrics and validated floating-point output bounds remain open. A customer evaluation must still establish that the supplied physical or operational model is appropriate and that total workflow costs improve on its actual repeated workload.


# DynaComp 0.6: sharp switching errors through a count lattice

Version 0.5 gave sharp box-error acceptance at a fixed operator and exact common quotients for uncertain families. Version 0.6 supplies approximate acceptance for a model whose operator and input map vary at every step through convex combinations of named endpoints. Its strongest new construction combines multi-affine extremal reduction, commuting-word compression and signed-support dynamic programming. The result is sharp for the declared switching/box domain. This is new machinery in this repository, not a claim that its constituent mathematical ideas are unprecedented.

## Why endpoints cover interior time-varying parameters

Suppose step j uses A[j] = sum alpha_s[j] A_s and B[j] = sum alpha_s[j] B_s, with nonnegative weights summing to one. A common projection uses G[j] and H[j] with those same weights. Hold the initial state, inputs and every other time's weights fixed. The full-minus-reduced output at a specified future time is affine in the weights at j: that time contributes its transition/input map once. Repeating this argument across all times shows multi-affinity in the independently chosen mode weights.

The maximum of an absolute affine function over a simplex occurs at a vertex. Applying this one time at a time shows that a worst case over all switching weights has a representative using endpoint modes only. Thus sharp evaluation of endpoint words covers every interior interpolation, including an arbitrary conductance value at every step. This does not justify checking only the two constant endpoint trajectories. Different endpoint orders can change the worst result. Nor does it turn the switching supremum into the sharp supremum for a parameter constrained to remain constant: that smaller domain can have a different extremal structure.

At a fixed endpoint word, the output error is linear in the initial state and all inputs. The componentwise box supremum follows from absolute coefficient sums, just as in 0.5. Initial and input deviations can be chosen independently to reinforce either sign of the central error. There is no need to enumerate their box corners in the algorithm. Exhaustive corner enumeration remains an independent test route on small instances.

## Compressing mode words without losing forcing order

For a suffix of a mode word, retain the paired row blocks a = C A_suffix and b = D G_suffix. When every pair of full operators commutes and every pair of projected operators commutes, these blocks depend only on the count of each mode in the suffix. A count vector becomes a node in a lattice, and adding one mode creates an edge. Both commutation conditions are checked exactly; projection can destroy commutation even when the original operators commute.

This compression alone is insufficient. The input error at a suffix node when the preceding mode is s is J = a B_s minus b H_s. Its contribution depends on which input time the edge represents and on the mode chosen there. Two words with the same final mode counts can therefore have different forcing contributions. A retained regression makes this explicit: two commuting shear operators with mode-dependent actuators give errors zero and one for the two orders of the same two modes. The dynamic program preserves this ordering information in its edge weights and predecessor path.

For each output, time t and orientation sigma in {+1,−1}, start at the zero-count node with score zero. At suffix depth d, the next edge corresponds to input time j = t−1−d. Give that edge weight sigma times J applied to the input center, plus the absolute coefficients of J applied to the input radii. Among paths arriving at a given count vector, retain the largest score and its predecessor word. Since future edge weights depend only on the current count and depth, discarding smaller scores is sound.

At a terminal count node, K = a minus b E is the initial-state error map. Add sigma times K applied to the initial center and its absolute coefficients applied to the initial radii. Maximize over terminal nodes and orientations. This yields the sharp componentwise error supremum for that output/time. Maximizing over the finite output/time set gives the sharp worst-horizon error. A stored predecessor path supplies the attaining endpoint mode word; coefficient signs supply the attaining initial field and heater/input word. Replay recalculates the bound and simulates that experiment to check its signed error.

For M modes, the number of count nodes through horizon T is the sum over d = 0, ..., T of binomial(d+M-1,M-1), equal to binomial(T+M,M). For two modes this is (T+1)(T+2)/2, versus 2^T endpoint words at the horizon. The number grows polynomially in T for fixed M, but can still be large with many modes. Evaluating every output time also repeats dynamic-program layers, and rational numerator/denominator growth remains a cost. Count-node reduction is not a constant-time or fixed-bit-complexity claim.

## Noncommuting cases and exact common closure

When either full or projected commutation fails, the implementation uses complete endpoint-word enumeration under an explicit symbolic-state budget. Exhausting that budget raises a resource-limit error rather than emitting an accepted partial maximum. The exact route remains practical for small horizons; noncommuting long-horizon sharp optimization needs further machinery.

An exact common quotient bypasses both routes. If DE = C, E A_s = G_s E and E B_s = H_s for every mode, every output error vanishes for every word. The certificate records zero bounds with one symbolic state. A regression exercises this on a noncommuting two-state family through eighty steps with a state budget of one. It represents 2^80 endpoint words without enumerating them because the identities establish all-future equality directly.

`switching_enrichment` starts from the target row space and optional seed directions. It forms the rational orthogonal lift V = Eᵀ(EEᵀ)⁻¹, projects every endpoint, and evaluates the sharp switching certificate. If tolerance fails and rank remains available, it adds an independent missing endpoint pullback. This is a general direction, not necessarily an original coordinate. Enrichment can reach common exact closure when the rank budget permits, but acceptance costs can hit the symbolic-state limit first. Its rank choices are heuristic, and error need not fall monotonically. It does not claim globally minimal approximate rank or minimum runtime.

## A smaller thermal model under changing conductance

The default switching example uses the 144-state 12-by-12 plate, corner-mean output, initial excess temperatures independently in [18,22], heater amounts independently in [0,2], and twenty steps. Conductance may vary independently at each step anywhere in [1/16,1/8]. These endpoint operators are affine functions of the same discrete Laplacian. Their orthogonal projections have the same affine structure, so both sets commute.

Nine retained directions give a sharp worst error of approximately 0.0828409143, below tolerance 0.1 across this entire changing-conductance domain. The calculation uses 231 count states while representing all 1,048,576 endpoint words at the horizon. In this case the worst bound equals the previous fixed-high-conductance bound, but the new calculation establishes that switching and interior values cannot make it larger. Equality of those two numbers alone would not prove that claim.

The separate common exact family still has seventeen states. At that rank the zero-error identity shortcut accepts arbitrarily long mode words without building a count lattice. Both results concern the supplied discrete model. They exclude continuum discretization error, spatially varying conductance, uncertain topology, physical installation mismatch and binary64 rounding.

## Reuse that survives input validation

`PreparedSwitching` is an optional numerical executor for two endpoint modes. The reduced path deep-copies the supplied plan, independently replays its rational certificate once, and prepares a numerically orthogonal QR basis and reduced endpoint operators. Subsequent batches reuse those operators. The full reference path validates the source and domain but does not incur the irrelevant reduction-certificate replay cost. Statistics distinguish preparation from batch execution.

A batch may share one initial field across many control plans. That field is encoded once, then its small reduced state is copied across scenarios. Per-scenario initial fields remain supported. Inputs have shape [steps, channels, scenarios], and interpolation parameters have shape [steps, scenarios] or a shared time sequence. Every parameter must be finite and between zero and one. The returned trace has shape [steps+1, outputs, scenarios]. Horizon and state/input box checks remain active by default.

Exact membership of binary64 data in rational boxes can be checked without constructing a Fraction for every quantity. For a rational lower endpoint, precompute the smallest binary64 value not below it. For an upper endpoint, precompute the largest binary64 value not above it. Comparisons against these inward-directed endpoints are equivalent to exact rational membership for finite binary64 values. Each endpoint is compared as a rational once, and `nextafter` corrects any outward rounding. The implementation handles subnormal and out-of-range endpoints, and a regression compares the result with exact Fraction tests around rational boundaries and extreme magnitudes.

This is an exact membership test for the input data, not interval-certified evolution. QR, interpolated matrix products and output computation still round. The rational bound concerns the ideal projected model. The application measures observed numerical agreement separately and does not fold a numerical rounding allowance into the claimed sharp algebraic supremum.

## Performance determines which representation is useful

The retained benchmark uses the production executor, domain checks enabled, a shared initial field and per-scenario time-varying interpolation. It compares full sparse evolution with a nine-state robust approximation over twenty steps, a seventeen-state exact quotient at the same horizon, and that exact quotient over four hundred steps. Every path returns every-time outputs. Five warmed measurements alternate execution order and retain raw timings. Generation, JSON I/O and installation are excluded; construction, source validation, certificate replay and numerical preparation are reported separately.

For 256 twenty-step scenarios, the approximate path takes about 0.000483 seconds versus full evolution at 0.00472 seconds, about 9.8 times faster in execution. Sharp robust construction costs about 8.70 seconds and reduced preparation/replay about 1.40 seconds, so estimated amortization requires about 2,382 batches. The exact seventeen-state path at that same horizon executes in about 0.000723 seconds, with about 0.558 seconds construction and 0.362 seconds reduced preparation. Its estimate is about 205 batches. The two short cases were measured in separate warmed trials; small timing variation should not be interpreted as a mathematical property.

For 256 four-hundred-step scenarios, full evolution takes about 0.0954 seconds and the exact reduced path about 0.0137 seconds, approximately seven times faster including domain checks. Construction costs about 0.528 seconds and reduced preparation about 0.356 seconds. Compared with full preparation and execution, the estimated break-even point is about 10.4 equally sized batches. The largest observed output difference is about 6.04 times 10 to the minus 14. All measured first compiled batches are slower than their full baselines after construction is included.

These results show why a smaller representation is not automatically the better deployment choice. The sharp approximate analysis is valuable for a tight state budget and changing-parameter assurance, but the common exact quotient can be substantially cheaper to build and already small enough to run efficiently. Timing medians are measurements, not guarantees about future runtime or customer profit. An actual workload needs its own initial-field sharing pattern, horizon, query count, integration costs and tolerance.

## Evidence and remaining fronts

The suite now has 131 methods. New tests compare the count-lattice result with every mode word and every state/input box corner on twelve random commuting models, check a mode-order forcing counterexample, compare noncommuting enumeration with direct word evaluation, exercise interior convex weights, verify the eighty-step exact shortcut, reject forged packets and incomplete budget results, test projection-induced noncommutation, check reusable/shared-field execution and directed membership, and round-trip the new CLI specification. All earlier retained receipts regenerate unchanged except additive summary entries.

The remaining mathematical fronts include noncommuting long-horizon sharp search, correlated input/state domains, approximate guarantees for nonlinear or spatially uncertain coefficient families, multivariate rank stratification and validated numerical evolution. The remaining product front is measured customer fit. This release makes changing-parameter assurance and reusable scenario execution concrete without claiming a universal speedup or an optimal representation search.


---

# DynaComp 0.7: all-time guarantees, exact stability and correlated domains

Versions 0.5 and 0.6 compute sharp worst-case reduction errors over boxes (and, in 0.6, under switching) for a finite horizon. That raises an obvious question it cannot answer: **does the guarantee survive after the horizon?** Version 0.7 answers it exactly for persistent input boxes. It supplies the stability machinery that answer needs, and it extends sharp bounds from boxes to correlated domains. As before, discovery may use floating point, but acceptance never does.

## The headline result

The 0.5 heated-plate workbench certifies a nine-state model of the 144-state plate. Its corner-mean error stays below 0.1 for twenty steps, for every initial field in [18,22] and every heater sequence in [0,2]. `examples/thermal_forever.py` keeps the same boxes and removes the horizon:

| Model | Twenty-step sharp bound | All-time worst error |
|---|---|---|
| Nine-state, finite-horizon enrichment (0.5) | 0.0828 | **refuted at t = 21** (error ≥ 0.105); bracketed in **[3.26665003, 3.26665003]**, attained near t = 121 |
| Seven-state, all-time enrichment (0.7) | — | **certified ≤ 0.0998 for every t ≥ 0** (bracket [0.0739, 0.0998]) |

The finite-horizon model was never wrong about its first twenty steps; it simply says nothing about later ones. Its error under a constant nominal heater converges to an exactly computed 0.389. With an adversarial bang-bang heater word, the error reaches 3.27. An independent binary64 simulation of the 144-state plate under the sign-pattern witness reproduces 3.26665003 to nine digits. The all-time enricher mixes target pullbacks with steady-state *moment* directions C(I − A)⁻ʲ, and with seven states it certifies the tolerance for every future time.

## All-time brackets

`all_time_bound` works on the exact observable quotient (E_o, G_o, H_o, D_o) of the full model, so the full and quotient outputs agree exactly. It stacks w = (E_o x, z) with A_w = diag(G_o, G), B_w = [H_o; H] and C_w = [D_o, −D]; the error is then e_t = C_w w_t. For a constant input centre a, the steady state is w* = (I − A_w)⁻¹B_w a and the steady error is e* = C_w w*. Let d₀ = lift(c) − w* and R_t = C_w A_wᵗ. The exact supremum at time t, over the initial box (radius r) and every input word in the box (radius b), is

S_t = |e* + R_t d₀| + Σ_k |R_t lift_k| r_k + Σ_{s<t} Σ_k |R_s B_w,k| b_k,

attained by a sign-pattern scenario. This is the 0.5 kernel identity, rewritten in the quotient.

**Sharp part (t ≤ T).** The rows R_t are propagated in dyadic arithmetic with 2⁻⁶⁴ resolution. Naive interval arithmetic would grow its radii like |A_w|ᵗ, which can diverge even for a stable A_w. DynaComp instead tracks the rounding error in the *dual* norm of a certified quadratic metric. If A_wᵀPA_w ⪯ γP, then A_w P⁻¹ A_wᵀ ⪯ γP⁻¹, so row errors contract by ρ = √γ, and the accumulated rounding error stays below its per-step size divided by (1 − ρ). Each functional R_t v is therefore enclosed within the computed value ± η_t‖v‖_P, and every term of S_t gets a rigorous upper and lower bound.

**Tail (t > T).** Every remaining term is bounded using |R_t v| ≤ κρᵗ‖v‖_P, where κ² = c P⁻¹ cᵀ. This gives

S_t ≤ |e*| + Σ_{s<T} |R_s B_w| b + κρᵀ (‖d₀‖_P + ν_x + β/(1 − ρ)),

which holds uniformly for all t ≥ T. The non-sharp part decays geometrically, so the bracket closes as T grows. The loop stops when it certifies the tolerance, when it refutes it, or when the non-sharp part falls below a slack (1e−6 by default).

**Refutation.** The lower bound is an attained value: it bounds S_t from below, and S_t is reached by an admissible scenario. If it exceeds the tolerance, the tolerance is false, not merely uncertified. The result is a statement about the model, independent of any proof technique.

**Repair.** `moment_rows` computes C(I − A)⁻ʲ on the exact quotient. Including the first moment in the encoder makes the reduced steady-state gain exact, so e* = 0. This follows from the standard left-Krylov moment-matching argument with any V such that EV = I, and a regression test checks it. `all_time_enrichment` lets the next moment and the next pullback compete each round, keeping whichever has the smaller lower bound. Candidates that fail are refuted early, so they are cheap to reject.

**Scope.** The bracket covers every t ≥ 0, every initial state in the box, and every input sequence that stays in a constant box forever. Both the quotient and the reduced dynamics must admit a certified metric; otherwise the result is `no_metric`, and finite-horizon certificates still apply. Replay checks P ≻ 0 and γP − A_wᵀPA_w ≻ 0 exactly, then recomputes the deterministic bracket. The floating-point Smith proposal is deterministic on a given platform. On another platform it may round differently, which would yield a different but equally valid P.

## Exact stability decisions and quadratic metrics

`stability.py` replaces the diagonal-only search of 0.1 with complete decisions and general metrics.

* **Discrete time.** `schur_stable` computes the characteristic polynomial over Q by Faddeev–LeVerrier and applies the Schur–Cohn recursion. The recursion reduces the polynomial to a primitive integer polynomial at each step, which keeps coefficients bounded. No eigenvalue is ever approximated.
* **Continuous time.** `hurwitz_stable` applies the Cayley transform (I − A)⁻¹(I + A), which maps the open left half-plane onto the open unit disk, and then the same recursion.
* **Certificates.** `quadratic_certificate` proposes P by floating-point Smith iteration on A/s. The scale s is chosen just above the estimated spectral radius, so that AᵀPA ⪯ s²(P − I) and γ approaches ρ(A)². The proposal is rationalized and accepted only through exact Sylvester (LDLᵀ pivot) tests of P ≻ 0 and γP − AᵀPA ≻ 0, with the smallest dyadic γ found by bisection. For the nonnormal matrix [[1/2, 8], [0, 1/2]], γ = 0.275 against an optimum of 0.25; the 0.1 diagonal search reported a factor above one.
* **Continuous certificates.** `continuous_certificate` checks AᵀP + PA ≺ 0 directly.
* **Switching.** `common_quadratic_certificate` searches for one P that contracts every mode, which proves stability under arbitrary switching. A pair of individually stable shear modes whose product is unstable is correctly left uncertified. Failure never claims instability.

Two independent routes check these claims. The randomized tests solve the Lyapunov equation exactly as an n²-variable rational linear system, and confirm the Lyapunov theorem's verdict (a unique positive definite solution exists exactly when the matrix is stable) agrees with Schur–Cohn and Routh-equivalent decisions on 280 random matrices. Every stable matrix in the certificate corpus received a certificate.

## Correlated domains

`sharp_zonotope_bound` replaces independent boxes with zonotopes over the whole scenario (x₀, u₀, …, u_{T−1}): a centre plus Σ ξ_g w_g with |ξ_g| ≤ 1. A generator may couple any coordinates and any times. Examples include one ambient offset shared by every cell, a heater held at one unknown level, or two heaters sharing a power budget. Because the error is linear in the scenario, its exact supremum is |L(centre)| + Σ_g |L(w_g)|, and the sign pattern that attains it is stored and executed during replay. Boxes are the special case of one generator per coordinate and time; a randomized test checks that the two give identical bounds, and another checks the result against exhaustive vertex enumeration.

Correlation should be reported honestly: in these diffusive examples the gain is modest. On the plate's twenty-step model, replacing independent ±2 cell uncertainty with a shared ±2 offset lowers the initial-state contribution from 0.0103 to 0.0072. The long-run failure above is driven by the mean initial field and by steady-state gain, not by input correlation, which is why moment directions repair it.

## Costs and limits

On the development machine, the thermal refutation took about 14 s, the converged bracket about 17 s (2,456 steps), all-time enrichment about 38 s, and independent replay about 18 s. Most of the time goes to the exact metric of the 17-state quotient, which is cached per model. The bracket needs constant input boxes; time-varying or correlated persistent domains, parameter-uncertain all-time bounds and continuous-time all-time bounds remain open. A common quadratic metric is sufficient for switching stability but not necessary. Nothing here validates the synthetic plate against a physical one.


---

# External engineering-model benchmark: DynaComp 0.7.1

The external benchmark uses published SLICOT system matrices recovered from the model-reduction source identified in our earlier work. The source is pinned to commit 98c1747030c0b1bf48d9353b2dcc9c4583867ca1. The selected models are a 200-state thin-rod heat equation, an 84-state PDE and a 48-state building-motion model. These are external engineering models, with controlled generated initial conditions and input sequences. They are not new measurements from operating equipment. Source identities, acquisition commands and conventions are in `benchmarks/README.md` and `benchmarks/slicot_manifest.json`.

Each continuous model is sampled through the augmented zero-order-hold matrix exponential. Every sampled binary64 coefficient is interpreted as an exact dyadic rational for certification; no decimal coefficient truncation is performed. The normalized output matrix has maximum row L1 norm one. The common contract is worst output error at most 0.01 over twenty discrete steps, for all independent initial-state and input coordinates in [-1,1]. The sample intervals differ by model and are retained explicitly. This contract concerns the sampled matrices, not continuous-time interpolation, physical calibration or matrix-exponential rounding.

The retained discovery method starts with the target and adds numerically proposed output-visible directions. Added directions are quantized to a dyadic denominator 2^20. Exact rational projection and sharp support evaluation decide whether each resulting four- or eight-state subspace meets tolerance. Independent replay reconstructs the projection, recomputes its box bound and executes the attaining corner. Numerical QR prepares execution in the same subspace; its roundoff is separate from the exact reduction claim. All six certificates replayed successfully. Three candidates meet tolerance and three do not.

For **heat-cont**, the four-state sharp bound is 0.0905726808 and the eight-state sharp bound is 0.000308047513. Their acceptance statuses are rank_budget_exhausted and tolerance_met. The eight-state repeated execution is 42.79 times faster than the full sampled model on this machine; that speed alone does not establish a useful result when acceptance fails. Its maximum observed error is 0.000308047513, with the retained attaining corner reproducing the rational supremum within numerical roundoff. The twelve-state balanced competitor has a binary64 box diagnostic of 3.11550533e-05 and execution speedup 33.23; that diagnostic is explicitly not an exact certificate.

For **pde**, the four-state sharp bound is 0.00270478545 and the eight-state sharp bound is 6.50011704e-07. Their acceptance statuses are tolerance_met and tolerance_met. The eight-state repeated execution is 12.10 times faster than the full sampled model on this machine; that speed alone does not establish a useful result when acceptance fails. Its maximum observed error is 6.50011704e-07, with the retained attaining corner reproducing the rational supremum within numerical roundoff. The twelve-state balanced competitor has a binary64 box diagnostic of 1.38727911e-08 and execution speedup 9.72; that diagnostic is explicitly not an exact certificate.

For **build**, the four-state sharp bound is 64.5600649 and the eight-state sharp bound is 12371.6251. Their acceptance statuses are rank_budget_exhausted and rank_budget_exhausted. The eight-state repeated execution is 5.24 times faster than the full sampled model on this machine; that speed alone does not establish a useful result when acceptance fails. Its maximum observed error is 12371.6251, with the retained attaining corner reproducing the rational supremum within numerical roundoff. The twelve-state balanced competitor has a binary64 box diagnostic of 43.707319 and execution speedup 4.06; that diagnostic is explicitly not an exact certificate.

The clearest positive case is the heat model: 200 states become eight while the rigorous normalized error stays below 0.000309, well inside the 0.01 contract. The repeated batch contains 256 scenarios and produces every output at every step. Full evolution takes 0.0101779 seconds per batch; the eight-state executor takes 0.000237878 seconds. Construction takes 32.816 seconds and independent replay takes 40.225 seconds. Including both, amortization requires approximately 7348 batches, or 1881140 twenty-step scenarios. A first batch is therefore dramatically slower with certification included. This external result replaces any suggestion that the synthetic plate's eleven-batch amortization should be expected generally.

The PDE model supplies another accepted case: four states achieve an exact normalized bound below 0.002705 and repeated execution is 15.27 times faster. Its construction plus replay costs 7.779 seconds and amortizes after about 3484 batches. Eight states improve its bound further but do not automatically improve economic value. A deployment should select a model according to the requested output tolerance, expected reuse and startup budget, rather than minimizing dimension alone.

The building-motion case is a substantive failure. Its four-state proposal has worst error about 64.56; eight states increase that to about 12,371.6. The attaining corners execute that failure in the full and reduced models. Adding directions is not monotonic, and an orthogonal Galerkin projection of a nonnormal stable system can have poor dynamics. These proposals are rejected, and their apparent runtime savings must not be sold as an accepted acceleration. The result does not prove that every four- or eight-state reduction is impossible. It identifies a weakness of the retained discovery heuristic. Stability-aware Petrov-Galerkin proposals and the new exact stability machinery are relevant next investigations, with no success claimed here.

The balanced-truncation comparison is an essential constraint on interpretation. Four-, eight- and twelve-state square-root balanced models are formed using discrete Lyapunov Gramians. Their numerical left-inverse residuals and Hankel singular values are retained. On the heat model, eight-state balanced truncation executes slightly faster than DynaComp but has observed initial-state-corner error exceeding 0.01; twelve states bring its numerical box diagnostic well below the threshold. Its zero-initial-state persistent-input response is substantially more accurate than DynaComp at matching ranks. This is consistent with balanced truncation optimizing an input-output objective rather than our arbitrary initial-state box objective. On the PDE model, the standard approach is also competitive. The experiment demonstrates accepted DynaComp reductions and a different accuracy contract; it does not establish broad superiority over conventional reduction.

Five warmed alternating timing trials use one BLAS thread for every engine and identical inputs. Random box scenarios are supplemented by a constant unit input from rest, a coherent unit initial field and both DynaComp worst-case corners. Per-batch data/domain checks, reduced matrix assembly, initial encoding and every-time output materialization are included. Common loading and sampling, network acquisition, JSON and gzip serialization are excluded from repeated execution. Sampling is reported separately; serialization and transfer further increase cold-start cost. A final timing sweep added the twelve-state balanced competitor and retained the original independent construction and replay measurements. Small timing differences are machine noise, not a mathematical property.

The exact cost is now visible rather than hypothetical. Rational coefficients in these dense sampled models can exceed Python's default decimal-string conversion cap; the trusted-data benchmark command enables their formatting locally. Six compressed packets total about 1.8 MB, although uncompressed rational text is much larger. Exact box acceptance adds tens of seconds to the heat model while ordinary numerical balanced discovery takes only a few hundredths of a second. Scalable rigorously enclosed dyadic kernels, reuse of prepared pullbacks and stability-aware candidate proposals are therefore more urgent than claiming additional warm-loop speed. They must be validated against these pinned models and adversarial corners.

All 148 tests pass after integrating the incoming v0.7 work. The six external packets independently replay, their attaining witnesses reproduce their sharp bounds numerically, and all observed DynaComp errors remain within the exact reduction bounds plus the explicit numerical comparison allowance. External data acquisition and expensive dense certification are intentionally outside routine CI. The complete release retains hashes, raw timings, rejected models and acquisition instructions, so a later run can report a genuine change rather than replace an inconvenient baseline.


---

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


---

# Infinite spectral sources and finite dynamics (DynaComp 0.9)

## What changes

Earlier DynaComp releases compare an approximate reduced model with a supplied finite model. That remains useful, but it leaves a preceding approximation unaccounted for: a differential equation may have infinitely many spatial modes, while the supplied matrix retains only a finite subset. Version 0.9 introduces explicit analytic source families and charges both omissions to the output error. It constructs a rational finite midpoint model, bounds its discrepancy from the infinite source, and can combine that bound with the existing acceptance calculation for any supplied finite reduction. A separate modal-selection route covers every future sample.

The sources are mathematical formulas. They are not field measurements or externally acquired engineering models. The heat source is a dimensionless one-dimensional diffusion equation, not the SLICOT 200-state heat benchmark. The Rogers–Ramanujan source uses the identities and coefficient bounds stated in K. Srinivasa Raghava, *Logarithmic moments of the Rogers–Ramanujan continued fraction and odd zeta values*, arXiv:2610.10578v1, https://arxiv.org/html/2610.10578v1. Those analytic identities are inputs to this implementation. Executable replay reconstructs coefficients and inequalities; it does not formally prove the paper or establish its results in Lean.

The central accounting is

\[
 |y_\infty-y_r|\le |y_\infty-y_K|+|y_K-y_r|.
\]

Here the first term includes omitted infinitely many modes and retained-coefficient enclosures. The second term concerns the actual rational matrices of the supplied smaller model. Neither numerical discovery nor a floating implementation is allowed to replace those matrices during checking. A source bound holding forever does not turn a finite-horizon second-stage bound into a permanent guarantee.

## Normalized Ramanujan moments as dynamics

The paper defines a normalized logarithmic moment \(B_G(s)\). Its Theorem 4 gives, for real \(s>1/2\),

\[
 B_G(s)=\sum_{N\ge1}c_NN^{-2s},\qquad
 |c_N|\le C(1+\log N),\qquad C=(\sqrt5+3)/2.
\]

For a specified starting integer \(n_0\ge1\), define one scalar state for every positive integer \(N\):

\[
 x_N(0)=N^{-2n_0},\qquad x_N(t+1)=N^{-2}x_N(t),
 \qquad y(t)=\sum_{N\ge1}c_Nx_N(t).
\]

Then \(y(t)=B_G(n_0+t)\). The normalized moment index is the discrete evolution variable; it is not physical time. This source is unforced. Passing a nonzero input to its executable model is rejected. The coefficient field is \(\mathbb Q(\sqrt5)\), and coefficient conversion is explicitly enclosed rather than silently rounded.

Equation (41) of the paper yields \(c_N=a_N+b_N\sqrt5\), with

\[
 b_N=\frac12\sum_{e\mid N}\frac{\chi(e)}e,
\]

\[
 a_N=-\frac12\sum_{e\mid N}\frac1e
 -\frac{\mathbf1_{5\mid N}}2\sum_{e\mid N/5}\frac{\chi(N/(5e))}{e}
 +\frac{\mathbf1_{25\mid N}}2\sum_{e\mid N/25}\frac1e.
\]

The implementation uses divisor sieves and exact fractions. The quadratic coefficient is strictly positive by its multiplicative local factors, so none of these coefficients vanishes. Tests independently expand the Rogers–Ramanujan Euler product and take the formal logarithm in the quadratic field; the first fifty coefficients agree exactly with the sieve calculation.

Let \(p=2n_0\). The function \((1+\log x)x^{-p}\) is decreasing for \(x\ge1\). Integral comparison therefore gives

\[
 \sum_{N>K}|c_N|N^{-2n_0}
 \le C K^{1-p}\left(\frac{1+\log K}{p-1}
                         +\frac1{(p-1)^2}\right).
\]

At an additional step \(t\), multiply this omitted-tail bound by \((K+1)^{-2t}\). Thus the supplied initial bound covers every future index and can become much smaller at later indices. The constants \(\sqrt5\) and \(\log K\) are enclosed by rational endpoints. For retained coefficients, midpoint/radius enclosures contribute the extra budget \(\sum_{N\le K}r_NN^{-2n_0}\), valid for every later index. Upward rounding of accumulated budgets keeps their representations small and preserves the inequalities.

Starting at \(n_0=2\), thirty-two modes give a source-only all-time upper bound approximately **0.000127809**. Selecting just the first mode and charging every removed retained contribution gives a one-state model with all-time upper bound **0.0140867**. Its output is constant near \(\phi^{-1}\). This is an approximation to the normalized sequence starting at order two; it is not an evaluation of every positive real moment, an approximation to unnormalized Gamma-growing moments, or an irrationality proof for odd zeta values. Starting at order one has a much less favorable conservative tail, and the source records that distinction.

## A general rational power-moment family

The same compiler also accepts the explicitly defined source

\[
 y(t)=\sum_{k\ge1}\frac{a}{k^{p+2(n_0+t)}},
 \qquad p\in\mathbb Z_{\ge0},\quad n_0\ge1,\quad a\in\mathbb Q.
\]

The retained coefficients and poles are already rational. The omitted initial contribution is bounded by

\[
 \frac{|a|}{(p+2n_0-1)K^{p+2n_0-1}}.
\]

This is a reusable family rather than a stored Ramanujan fixture. Signed amplitudes and the identically zero amplitude are handled explicitly. Additional infinite sources require their own formulas, convergence domains and derived tails. The interface does not accept an arbitrary user-asserted remainder as evidence.

## The infinite heat rod

Consider the dimensionless boundary-value problem

\[
 \partial_t\theta=\partial_{xx}\theta+u(t),\quad 0<x<\pi,
 \qquad\theta(0,t)=\theta(\pi,t)=0,
 \qquad\theta(x,0)=a.
\]

The target is the spatial mean \(y(t)=\pi^{-1}\int_0^\pi\theta(x,t)\,dx\). The constant initial field is understood in the usual integral/spectral sense; its endpoint mismatch is harmless for this target. Fourier sine expansion of the uniform initial field and uniform forcing gives only odd modes. With \(k=1,3,5,\ldots\),

\[
 \dot z_k=-k^2z_k+u(t),\quad z_k(0)=a,
 \qquad y(t)=\frac8{\pi^2}\sum_{k\text{ odd}}\frac{z_k(t)}{k^2}.
\]

The declared domains are \(|a|\le X\) and \(|u(t)|\le U\). The initial modes are correlated because they share the same scalar \(a\). The scalar estimate

\[
 |z_k(t)|\le X e^{-k^2t}+U/k^2
\]

holds for measurable bounded forcing. If the first \(m\) odd modes are retained, convexity on cells centered on the omitted odd integers gives

\[
 \sum_{j\ge m}(2j+1)^{-p}
 \le \frac{(2m)^{1-p}}{2(p-1)},\qquad p>1.
\]

Consequently the omitted output is bounded at every continuous time by

\[
 \frac8{\pi^2}\left[
   \frac{X e^{-(2m+1)^2t}}{4m}+\frac{U}{48m^3}\right].
\]

This is the bound for omitted modes, not a continuous-time guarantee for an arbitrary subsequently projected model. The finite executable model samples at rational interval \(\Delta>0\) and assumes the heater input is held constant over each interval. Its exact scalar coefficients are

\[
 a_k=e^{-k^2\Delta},\qquad b_k=(1-a_k)/k^2,
 \qquad c_k=8/(\pi^2 k^2).
\]

Both the exponential poles and the common output factor are enclosed. Rational midpoints define the actual stored finite model. Defining \(b_k\) from the midpoint pole preserves the exact scalar steady-state gain \(1/k^2\); rounding pole and input coefficient independently would lose that property.

## Retained coefficient and sampling errors

Let \(\hat a_k\) be a midpoint pole, \(\delta_k\) its enclosing radius, and \(\rho_k\) the upper pole endpoint. The compiler requires \(\rho_k<1\). If insufficient precision cannot separate a very slow mode from one, it raises an explicit error rather than divides by a vanishing decay margin. The true and midpoint filters have tied input coefficients \((1-a_k)/k^2\) and \((1-\hat a_k)/k^2\).

For \(0\le a_k,\hat a_k\le\rho_k<1\),

\[
 |a_k^j-\hat a_k^j|\le j\rho_k^{j-1}\delta_k.
\]

The initial-state discrepancy is therefore bounded for every sample by \(X\delta_k/(1-\rho_k)\). Splitting the forced impulse difference into its coefficient and power differences and summing geometric series yields

\[
 |z_k-\hat z_k|\le
 \frac{X\delta_k}{1-\rho_k}
 +\frac{U\delta_k}{k^2}
 \left[\frac1{1-\rho_k}
       +\frac{1-\hat a_k}{(1-\rho_k)^2}\right].
\]

If the output coefficient midpoint is \(\hat c_k\), with radius \(r_k\), the retained-mode output discrepancy is at most

\[
 r_k(X+U/k^2)+|\hat c_k|\,|z_k-\hat z_k|.
\]

Sum those terms and add the two omitted-mode contributions. This covers every future sample under every permitted held-input sequence. The assumptions are bounds on a prescribed mathematical rod; they do not validate conductivity, geometry or operating ranges against physical measurements.

## All-time selection and arbitrary second-stage reduction

`compress_spectral_modes` keeps any distinct nonempty list of modal coordinates, including a reordered list. For each deleted rational mode its output weight is multiplied by an all-time state-size bound. For unforced moments that size is its prescribed initial value. For the heat model it is \(X+U/k^2\). Adding the deleted contributions to the source budget gives an all-time bound for the selected model. This route does not need an exact observable quotient or a dense Lyapunov solve.

A general projection is handled by `compose_spectral`. Its actual encoder, transition, input and output matrices go through the existing outward-dyadic box calculation. For heat, the common initial field is enlarged to the independent modal box \([-X,X]^m\); this loses some sharpness but keeps acceptance valid. The supplied input samples have independent intervals \([-U,U]\). The source bound and finite-reduction upper bound are added. A failed lower witness for the enclosing finite model does not automatically refute the correlated infinite source, so the overall unsuccessful status is **not_certified**, not **refuted**.

The retained heat demonstration uses \(\Delta=0.1\), \(U=1\), and tolerance **0.02**. With \(X=1\), sixty-four midpoint modes give source-only upper bound **0.00316635**. Selecting eleven modes gives all-time upper bound **0.0184221**. Optional Lyapunov-coordinate discovery proposes one, two and three states; their combined twenty-sample upper bounds are approximately **0.193024**, **0.0381488** and **0.00466571**. Only the three-state candidate meets the tolerance. Its reduced Euclidean contraction is checked exactly, but its output-error acceptance still stops at sample twenty.

For the different domain \(X=0\), selecting one mode from thirty-two gives an all-time upper bound **0.011897564537**. A unit constant input approaches an actual error

\[
 \pi^2/12-8/\pi^2\approx0.011897564285.
\]

An independent Decimal calculation reproduces this limit, showing that the conservative receipt is within about \(2.52\times10^{-10}\) in this particular case. This near-sharpness concerns the one-mode, from-rest rod and should not be generalized to arbitrary selected models. The two initial domains must not be pooled into one performance claim.

## Two-cusp evaluation

`ramanujan_value(t)` encloses \(R(e^{-t})\) for positive rational \(t\). At \(t\ge6\), it uses the direct Euler product. At smaller arguments, it applies the classical reciprocal identity and evaluates at \(t'=4\pi^2/t\). The threshold is an implementation choice, not a singularity or asymptotic extrapolation.

For a direct product truncated after \(J\), with \(q=e^{-t}\), the omitted logarithm is bounded in absolute value by

\[
 \sum_{j>J}\frac{q^j}{1-q^j}
 \le \frac{q^{J+1}}{(1-q)(1-q^{J+1})}.
\]

The product endpoints are multiplied by the corresponding exponential factors. All constant evaluations and successive product operations are enclosed outward. Reciprocity is rearranged as

\[
 R(e^{-t})=\frac{1-\phi R(e^{-4\pi^2/t})}
                 {\phi+R(e^{-4\pi^2/t})},
\]

avoiding a subtraction of nearby large terms. At \(t=10^{-6}\), the reciprocal route needs just one product term and returns an interval about **3.31e-24** wide at eighty bits. This is a value enclosure, not an implementation of the ordinary two-cusp Mellin integrals. The paper's incomplete-gamma moment continuation remains a separate extension.

## Rational enclosures and execution

Square-root endpoints come from integer square roots on a dyadic grid. Positive rational logarithms use power-of-two range reduction and the atanh series with an explicit geometric remainder. Machin's identity expresses \(\pi\) through two alternating arctangent series. Negative exponentials use range reduction, alternating Taylor remainders and outward-rounded repeated squaring. Series stop when their proved remainders meet the internal budget. No approximate numerical special function supplies an accepted endpoint.

`simulate_spectral` executes the actual stored matrices in rational arithmetic. It checks the source input domain, the common initial rod field and any second-stage horizon. Uncertified packets are rejected by default; an explicit override allows diagnostic execution without changing their status. The mathematical bounds do not include runtime floating-point evolution error. The exact executor is suitable for checking and examples, not a high-throughput production claim.

Run `python examples/infinite_spectral.py` to regenerate receipts and the summary. The analytic sources and modal-selection routes use only the standard library. Installing the numerical extras enables the general-direction heat proposals. CLI kinds are `spectral`, `spectral_subset`, `spectral_composition` and `ramanujan_value`; their example JSON specifications all compile and replay. The budgeted geometric search returns its final packet and tested ranks; exhausting its mode budget never becomes acceptance and no smallest-rank claim is made.

## Exact limits and validation

The Ramanujan generating function has a nonzero pole at every positive square according to the cited theorem. A fixed finite-dimensional time-invariant linear sequence realization has a rational generating function with only finitely many poles. Thus exact realization of the full sequence by that class is impossible. The heat transfer function similarly has infinitely many visible odd modal poles. These observations motivate approximation with explicit tails; they do not exclude nonlinear evaluation algorithms, time-varying constructions or other forms of compact description.

The new tests check integer endpoint inequalities, high-precision Decimal constant values, independent continued-fraction evaluation, formal quadratic-field coefficient expansions, long independent tail sums, source-domain guards, exact executions, finite/infinite error composition, packet tampering and precision insufficiency. The heat transient test independently evolves three hundred exact-analytic scalar modes with eighty-digit Decimal arithmetic and charges its own omitted reference tail. The optional projection test executes both stored models exactly and verifies their discrepancy against the accepted finite bound. These are computational cross-checks in addition to the stated mathematical derivations.

Open directions are descriptor/mass-matrix sources, calibrated engineering workloads, correlated second-stage initial domains, general continuous-time projected error bounds, other analytic spectral families and efficient accepted floating execution. This release supplies explicit source-to-model guarantees for three defined families, not a universal compiler for arbitrary infinite-dimensional equations.


The complete release suite has **188 tests passing with numerical extras**. With only the standard library available, **170 tests pass and eighteen optional tests skip**. Twenty-nine new spectral/evaluation test methods extend the earlier 159-test suite. Every example specification also passes the compile, JSON round-trip and replay route. The analytic demonstration emits accepted and unaccepted candidates, preserving their distinct statuses rather than presenting all proposed reductions as successful.
