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
