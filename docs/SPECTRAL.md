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
