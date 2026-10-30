# From mathematical reduction to a measured executable workload

**DynaComp 0.11.0, 10 October 2026.** This extension supplies a reusable external-workload interface, a conventional compiled binary64 executor with an explicit error contract, and measurements against full dense BLAS and square-root balanced truncation. It changes the practical question from whether a small model can be constructed to whether an accepted executable model is useful after paying for its construction. The answer is positive on two retained external sources and negative on a third. All three outcomes are part of the release.

## The source and the operating question

The public SLICOT model-reduction collection includes a 200-state heat model, an 84-state convection/reaction model, and a 48-state building model. Their original MAT files and repository revision were pinned during the earlier external benchmark. This release automatically uses the retained sampled and normalized-output rational source in each packet and checks its SHA-256 against `benchmarks/slicot_manifest.json`. Upstream is `https://github.com/SLICOT/Benchmark-ModelReduction`, revision `98c1747030c0b1bf48d9353b2dcc9c4583867ca1`. The official collection description is at `https://www.slicot.org/20-site`; it describes continuous descriptor systems and documents the source dimensions and applications.

The mathematical source here is the stored sampled model, not the original differential equation. Every finite binary64 sampling coefficient was imported as its exact dyadic rational. Sampling error from the matrix exponential and agreement with a measured physical system are outside this guarantee. Output rows were normalized in the earlier benchmark, so the tolerance is expressed in those normalized units. The importer for new customer models preserves output units unless normalization is explicitly selected, and records its divisor.

The retained operating contract permits every independent initial coordinate in [-1,1], every input coordinate in [-1,1] at each of twenty transitions, and requires absolute output error at most 0.01 at time zero and after each transition. These are experimental operating boxes, not measured customer trajectories. A claim about this contract does not imply a claim about arbitrary continuous trajectories, arbitrary parameter uncertainty, sensor calibration, or a real building's plausible initial conditions. A supplied customer model should use a domain with a meaningful engineering interpretation.

## Separate products and additions

The v0.10 executor accumulates exact integer products of stored dyadic operands and rounds once per row. It is portable and useful for exact numerical semantics, but its Python integer machinery is an expensive path to a high-throughput scenario engine. This release adds a different protocol: each dot starts at positive zero, forms one binary64 product at a time, and adds it to the accumulator in left-to-right order. The state transition concatenates the reduced-state terms and input terms into the same dot. Encoding, transition and output each follow this rule.

The C implementation is shipped as source in the installed package. It is compiled once per process, using a GCC-compatible compiler, `-fno-fast-math`, `-ffp-contract=off`, and `-frounding-math`. Volatile product and accumulator stores make the separate arithmetic operations explicit. The implementation does not delegate to BLAS and does not permit fused multiply-add or reassociated summation. A compiler/preprocessor check requires binary radix, 53-bit precision, the binary64 exponent range, and directly evaluated double arithmetic. Execution checks round-to-nearest mode and gradual underflow before accepting a batch. The underflow checks compare stored integer encodings, because floating comparisons in a denormals-are-zero environment could otherwise mistake an expected subnormal for zero.

This is an explicit executable contract, not a certificate for arbitrary optimized C, NumPy, a GPU kernel, or a customer's existing BLAS library. It assumes the C compiler and platform implement the stated operations. No compiler correctness theorem, Lean proof of C semantics, or new kernel-checked bridge is claimed. A source hash and compiler version are recorded in the benchmark report. Changing the arithmetic engine requires a corresponding execution analysis.

## A dot bound including underflow

Let u = 2^-53 be binary64 unit roundoff and eta = 2^-1075 be half the least positive subnormal. Subject to finite rounding, one scalar operation obeys the mixed bound |fl(v)-v| <= u|v| + eta. The normalized rounding convention is described in the NIST DLMF, Section 3.1 (`https://dlmf.nist.gov/3.1`). The additional eta term is retained here so the analysis also covers gradual underflow instead of assuming every result is normal.

For a length-k dot, let M bound the sum of absolute exact products of its stored operands. There are at most 2k scalar rounding operations: k multiplications and k accumulator additions, including the initial addition to zero. Define L = 2k and require Lu < 1. Repeated substitution of the scalar rounding model bounds the coefficient perturbations by `(1+u)^L-1` and the propagated additive underflow terms by at most `L*eta*(1+u)^L`. The elementary inequality `(1+u)^L <= 1/(1-Lu)` therefore gives the conservative sufficient error budget

`gamma_L * M + L*eta/(1-Lu)`, where `gamma_L = Lu/(1-Lu)`.

The implementation evaluates this expression in exact rational arithmetic and rounds the positive bound upward to a compact dyadic value. A zero magnitude gives a zero arithmetic error. It also requires M plus the unrounded error expression to remain below the largest finite binary64 value, which conservatively excludes overflow of products and intermediate sums. Domains where cancellation would make the final result safe can still fail this absolute-magnitude condition. Failure means the supplied sufficient bound could not accept the domain; it is not a theorem that every execution overflows.

This bound is intentionally more conservative than a specialized sharp dot analysis. Its advantage is a simple operation count, an additive underflow treatment, and direct compatibility with the existing state-error recurrence. The randomized tests compare native results against a separate implementation that converts stored words to exact integers/rationals, rounds each product independently, and rounds each sum independently. This checks the intended sequence of operations, rather than comparing native C with another call to the same C routine.

## Loading, state propagation and composition

Let E,G,H,D be the supplied rational encoder, transition, forcing and output matrices. Hat matrices denote their loaded binary64 values, interpreted exactly. Let Z bound the ideal reduced state and e bound the difference between the actual stored native state and that ideal state. Coefficient-loading errors are retained separately. The initial error is bounded by `|Ehat-E|X + dot_budget(|Ehat|X,n)`, where X is the initial absolute state envelope.

At each transition, the ideal envelope advances with `Z_next = |G|Z + |H|U`. The actual deviation is bounded by `|Ghat|e + |Ghat-G|Z + |Hhat-H|U + dot_budget(|Ghat|(Z+e)+|Hhat|U,r+m)`. The output deviation is bounded by `|Dhat|e + |Dhat-D|Z + dot_budget(|Dhat|(Z+e),r)`. Every envelope is rounded outward; n is source dimension, r is reduced dimension, and m is input-channel count. The execution certificate records every output and state bound across the declared finite horizon.

The model-reduction certificate independently bounds output error between the exact rational source and ideal reduced equations. Sequential composition replays that source certificate and adds the execution error at the same output and time. Its maximum is a maximum of these time-matched sums, not necessarily the sum of the two unrelated global maxima. This can be tighter while remaining conservative. Acceptance requires the resulting combined upper bound to meet the target and any requested reduced contraction condition to hold.

Two new packet schemas, `dynacomp.sequential-binary64/1` and `dynacomp.sequential-composition/1`, reconstruct from their source matrices and domains. Existing integer-accumulated packet formats and receipts are unchanged. The general JSON CLI supports `sequential_execution` and `sequential_composition`; mathematical compilation and replay do not require NumPy or a compiler. Preparing and running the batched native arrays does require the numerical extras and a supported compiler.

## Searching two methods against one contract

The pilot interface proposes candidate reduced models in ascending supplied rank order. One method uses the existing Lyapunov-coordinate observable enrichment. The other is a transparent square-root balanced truncation implementation using SciPy discrete Lyapunov solves, eigenvector Gramian roots and an SVD. The latter is an established algorithmic baseline, but these timings are not a measurement of the SLICOT binary toolbox. Numerical Gramian failures or insufficient Hankel rank are reported rather than interpreted as exact mathematical rank facts.

Both candidates enter the same outward-dyadic finite-box acceptance procedure. An obviously failed model does not pay for a separate native composition. A model passing the ideal rational target must additionally pass the conventional execution bound. Both accepted reductions subsequently use the same C executor. The observable method also requires exact reduced Euclidean contraction; balanced acceptance in this experiment is finite-horizon and does not impose that extra norm condition. This distinction is retained in packets, rather than silently treating a finite-horizon guarantee as an all-time stability certificate.

The selected rank is the first passing rank along each method's own candidate path. It is not a proof of minimum dimension among all possible reduced models. The comparison matches a common maximum admissible error, 0.01, rather than forcing equal realized errors. In particular, the accepted heat balanced candidate has a substantially smaller actual bound than the accepted observable candidate. The ordinary balanced objective concerns zero-initial-state input/output behavior; arbitrary initial-coordinate boxes impose an additional demand. That helps explain why a balanced proposal may fail this experiment without implying that balanced truncation is generally inferior.

## Saved plans and real supplied scenarios

`dynacomp-pilot compile` accepts a rational discrete workload and emits a plan with both rank histories, accepted packet artifacts, hashes, source identity, operating boxes and setup measurements. `dynacomp-pilot evaluate` adds execution measurements and an interactive HTML report. `dynacomp-pilot run` re-evaluates a saved plan. `dynacomp-pilot simulate` accepts actual initial and input arrays from `.npy` files and writes the complete output trace. This provides a practical route for a prospective user to evaluate their own repeated workload rather than reading a research script.

Saved packet paths must remain inside the plan directory. Artifact hashes and summary fields are checked, and independently replayed packets must match the workload's source, tolerance, initial domain and input domain. Retained dense plan and workload JSON can be gzip-compressed without changing the interface. Packet properties return copies; internal loaded matrices are read-only. A batch can share one initial field or provide one per scenario. Stored binary64 inputs must lie inside inward-rounded exact endpoints, and nonfinite data or a different horizon is refused. Loading arrays uses `allow_pickle=False`.

`import-mat` performs automatic ingestion for identity-mass, zero-feedthrough real matrices. A supplied positive rational sample interval triggers zero-order-hold sampling through an augmented matrix exponential. Omitting the interval explicitly treats the input as discrete-time. General descriptor matrices, nonzero direct feedthrough, complex matrices and numerical sampling error certification are open; unsupported cases are refused. Default boxes are placeholders for exploration and must be changed when a user's real operating domain differs.

## Retained measurements and their limits

The release measurements use 256 scenarios per batch, twenty steps, seven warmed trials with alternating method order, and one BLAS thread. Every timed method includes domain checking, initial encoding where applicable, and materialization of the complete output array. The full reference is a prepared dense NumPy BLAS simulation. Sampled differences against that reference are diagnostic only: its floating evolution is not certified. The reduced acceptance guarantee is independently established against the rational sampled source and does not depend on trusting BLAS output as a proof oracle.

For the 200-state heat source, observable enrichment accepts seven states with combined upper bound approximately 0.006988834283. Its global execution-only bound is approximately 9.63e-13. The retained native execution is 5.89 times faster than full BLAS after setup. Balanced truncation accepts nine states with combined bound approximately 0.000599539489 and warm speedup 4.71. The conservative in-process setup totals are about 2.45 seconds and 2.87 seconds respectively. Their estimated break-even points are roughly 222 and 283 batches of 256 scenarios. At one thousand such batches, projected setup-inclusive speedups are 2.63 and 2.18, respectively. These totals extrapolate measured medians rather than recording a thousand-batch sustained run.

For the 84-state convection source, both methods accept four states. Observable and balanced combined bounds are approximately 0.002163652598 and 0.003059127706. Warm speedups are 5.99 and 6.06, so this is not evidence of a meaningful observable execution advantage over the balanced baseline. Estimated break-even points are roughly 212 and 237 batches; projected thousand-batch total speedups are 2.78 and 2.64. Differences this small are sensitive to timing variation and setup accounting.

For the 48-state building source, observable enrichment accepts thirteen states with combined bound approximately 0.000901406147. Conventional execution rounding is noticeably larger here: its global envelope is approximately 1.88e-7, still below the tolerance margin. The accepted native engine runs at only 0.45 times the full BLAS speed, so it is slower and has no measured break-even. Balanced ranks one through sixteen fail the supplied initial/input box target. This result demonstrates that dimension reduction and certified accuracy do not imply a runtime benefit. It also keeps the failed baseline visible instead of discarding the source from the report.

Setup accounting includes in-process source conversion, proposal search, certificate construction, independent replay, native compilation and preparation. Both standalone reduced methods are charged compilation even if the process already cached the library. Search shares source pullbacks; the first method pays their initial construction. For a conservative standalone estimate, balanced setup is charged an additional empty-cache replay to cover the shared source pass. This can overestimate its standalone cost and is reported explicitly. Interpreter startup, dependency imports, file I/O, upstream acquisition and matrix-exponential sampling are excluded. A deployment comparison should add those costs when they differ between routes.

## What remains commercially unproved

The interface and evidence are sufficient to run a concrete customer evaluation. They are not evidence of a paying customer, safety certification, universally superior reduction, or an industrial throughput win. The next empirical step is a customer's real source matrices, operating ranges, expected number of repeated batches and a baseline chosen with that customer. The product should refuse or recommend the full engine when setup cannot amortize or the smaller executor is slower. The retained building failure is an example of precisely that situation.

Continuous intersample guarantees, all-time conventional rounding invariants, correlated-domain native composition, descriptor models and an optimized conventional batch kernel remain open. Future optimized kernels must retain their operation contract or obtain a new rounding analysis. The important advance of this release is narrower and concrete: a supplied sampled workload can now travel through reduction, independently checkable model error, actual conventional execution error, guarded saved-plan execution and a transparent cost comparison.

## Windows native portability

The Windows path targets **x86-64 CPython and MinGW-w64 GCC (UCRT64)**.
MSVC, 32-bit Windows, ARM64 and MSYS/Cygwin compilers are not supported by this
path. Install the numerical extras with the Windows Python interpreter:

```powershell
python -m pip install -e ".[numerical]"
$env:DYNACOMP_CC = 'C:\msys64\ucrt64\bin\gcc.exe'
$env:PATH = 'C:\msys64\ucrt64\bin;' + $env:PATH
python -m unittest discover -s tests -p test_native.py -v
dynacomp-pilot evaluate applications/customer_pilot/two_state.json --output pilot-windows --queries 8 --trials 3
python tools/check_native_portability.py pilot-windows
```

The same compiler path can be set in Git Bash with
`export DYNACOMP_CC='C:/msys64/ucrt64/bin/gcc.exe'`. Use Windows Python, not an
MSYS Python. `DYNACOMP_CC` is one executable path, not a command with flags.
Absent this override, Windows discovery tries `x86_64-w64-mingw32-gcc` and
`gcc`, and verifies an x86-64 MinGW target. Linux retains cc/gcc/clang discovery
and its existing shared-library flags.

Windows builds a DLL, explicitly exports the cdecl `run` entry point, and loads
it with `ctypes.CDLL` by absolute path. SSE2 evaluation and static libgcc avoid
x87 excess precision and an extra libgcc runtime DLL dependency. Keep the UCRT64
bin directory on PATH for compiler tools. At process exit, the library is
unloaded before temporary build files are deleted. Builds and native execution
remain explicit; certificate construction/replay still needs no compiler.

Both platforms use C11, separate volatile products and sums, no fast math,
no FP contraction, and rounding-aware compilation. The preprocessor checks and
per-batch RN-even/subnormal guards remain enabled. The runtime refuses an
unsupported floating-point environment; it does not change the caller's mode.
No certificate schemas or error formulas change.

CI adds Windows Server 2022 x64 / CPython 3.11 / MSYS2 UCRT64 GCC and Ubuntu
24.04 / CPython 3.11, alongside the established Linux matrix. It requires
successful native compilation before running tests, so missing compiler support
cannot silently turn native coverage into skips. Target-compiled fenv constants
replace Linux-specific rounding constants in tests. Tests cover halfway ties,
subnormal generation/transport, FTZ/DAZ rejection, overflow refusal, and random
traces against independent integer/rational rounding.

Each runner evaluates the supplied two-state customer example and retains its
plan, packets, timings and a `portability.json` report. Fixed endpoint/interior
stored-input scenarios are independently evolved using exact rational **full
source** equations. Every native output must satisfy its time-matched combined
model/execution bound. Repeated executions on each runner must agree. A dependent
job compares the two platforms against the sum of their certified bounds;
platform-dependent reduction coefficients may differ. Hardware timings are never
required to match. These sampled checks supplement the existing all-box
certificates rather than replacing them.

A configured CI environment is not evidence of a passing environment. Consult
the PR's actual Actions jobs and retained compiler/build metadata before claiming
Windows validation; Windows 11 on a customer laptop requires its own run.
