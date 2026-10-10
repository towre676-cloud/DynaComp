"""Infinite spectral families -> rational finite dynamics with analytic tails.

Supported source semantics are explicit formulas, not arbitrary asserted tail
bounds. Ramanujan identities are mathematical inputs from the cited paper;
replay checks executable arithmetic and bounds, not a formal proof of that paper.
"""
from fractions import Fraction
from functools import lru_cache
from .rational import q, mv, matrix, addv
from .enclosures import sqrt_integer, log_rational, pi_interval, exp_negative, precision
from .uniform import SparseLinearModel
from .validated import validated_box, replay_validated

PAPER = "https://arxiv.org/html/2610.10578v1"


def _upper(value, bits):
    """Round a nonnegative budget upward to keep receipts bounded in size."""
    scale = 1 << bits
    return Fraction((q(value)*scale).__ceil__(), scale)


def _integer(value, name, low=1, high=4096):
    if type(value) is not int or not low <= value <= high:
        raise ValueError(f"{name} must be an integer in [{low},{high}]")
    return value


def _nonnegative(value):
    value = q(value)
    if value < 0:
        raise ValueError("Nonnegative amplitude/domain required")
    return value


def character5(n):
    return (0, 1, -1, -1, 1)[n % 5]


def ramanujan_coefficients(count):
    """Return exact (a,b) with c_N = a + b sqrt(5), N=1,...,count.

    Divisor sieves implement equation (41), independently of continued-fraction
    evaluation. Fractions retain all cancellations and coefficient signs.
    """
    _integer(count, "Coefficient count")
    harmonic = [Fraction(0) for _ in range(count+1)]
    direct = harmonic.copy()
    reverse = harmonic.copy()
    for divisor in range(1, count+1):
        inv = Fraction(1, divisor)
        for n in range(divisor, count+1, divisor):
            harmonic[n] += inv
            direct[n] += character5(divisor)*inv
            reverse[n] += character5(n//divisor)*inv
    return tuple((-harmonic[n]/2
                  - (reverse[n//5]/2 if n % 5 == 0 else 0)
                  + (harmonic[n//25]/2 if n % 25 == 0 else 0),
                  direct[n]/2) for n in range(1, count+1))


def _field_interval(a, b, root):
    values = (a+b*root[0], a+b*root[1])
    return min(values), max(values)


def _midpoint(interval):
    lo, hi = interval
    return (lo+hi)/2, (hi-lo)/2


@lru_cache(maxsize=32)
def _constants(bits):
    root = sqrt_integer(5, bits)
    pi = pi_interval(bits)
    return root, (8/pi[1]**2, 8/pi[0]**2)


def _source(family):
    """Canonicalize source formulas; reject unnoticed/unsupported parameters."""
    kind = family.get("kind")
    fields = {
        "ramanujan_log": {"kind", "start"},
        "power_moments": {"kind", "start", "power", "amplitude"},
        "heat_rod": {"kind", "dt", "initial_amplitude", "input_bound"},
    }
    if kind not in fields or set(family)-fields[kind]:
        raise ValueError("Unsupported spectral family or parameter")
    if kind == "heat_rod":
        dt = q(family.get("dt", "1/10"))
        if dt <= 0:
            raise ValueError("Positive sample interval required")
        return {"kind": kind, "dt": str(dt),
                "initial_amplitude": str(_nonnegative(family.get("initial_amplitude", 1))),
                "input_bound": str(_nonnegative(family.get("input_bound", 1)))}
    result = {"kind": kind, "start": _integer(family.get("start", 1), "Starting moment", high=128)}
    if kind == "power_moments":
        result.update(power=_integer(family.get("power", 0), "Coefficient power", low=0, high=128),
                      amplitude=str(q(family.get("amplitude", 1))))
    return result


def _build(family, modes, bits):
    """Build midpoint dynamics and independently accounted all-time budgets."""
    family = _source(family)
    _integer(modes, "Modes")
    precision(bits)
    kind = family["kind"]
    entries, b, c, center, radius, modal = [], [], [], [], [], []
    coefficient_error = Fraction(0)
    if kind in ("ramanujan_log", "power_moments"):
        start = family["start"]
        if kind == "ramanujan_log":
            root = sqrt_integer(5, bits)
            pairs = ramanujan_coefficients(modes)
            weights = [_field_interval(a, b, root) for a, b in pairs]
            exponent = 2*start
            # Integral comparison for a decreasing log-weighted p-series.
            log_hi = log_rational(modes, bits)[1]
            tail = (root[1]+3)/2*Fraction(1, modes**(exponent-1))*((1+log_hi)/(exponent-1)+Fraction(1, (exponent-1)**2))
            obstruction = "Assuming the cited Theorem 4: infinitely many nonzero poles; no exact fixed finite-dimensional time-invariant linear realization of the full moment sequence."
        else:
            amplitude, power = q(family["amplitude"]), family["power"]
            weights = [(amplitude/Fraction(k**power),)*2 for k in range(1, modes+1)]
            exponent = power+2*start
            tail = abs(amplitude)*Fraction(1, (exponent-1)*modes**(exponent-1))
            obstruction = "Nonzero amplitude gives infinitely many distinct poles; zero amplitude is identically zero."
        for i, weight in enumerate(weights):
            k = i+1
            pole = Fraction(1, k*k)
            mid, rad = _midpoint(weight)
            initial = pole**start
            entries.append([i, i, str(pole)])
            b.append(["0"])
            c.append(str(mid))
            center.append(str(initial))
            radius.append("0")
            coefficient_error = _upper(coefficient_error+rad*initial, bits)
            modal.append({"index": k, "pole_interval": [str(pole)]*2,
                          "weight_interval": [str(v) for v in weight]})
        domain = {"initial_center": center, "initial_radius": radius, "input_bound": "0"}
        budgets = {"omitted_initial": str(_upper(tail, bits)), "omitted_forcing": "0",
                   "retained_coefficient_and_sampling": str(coefficient_error)}
        semantics = "Unforced normalized moment at index start+t, every integer t>=0; input must be zero."
    else:
        dt = q(family["dt"])
        initial, forcing = q(family["initial_amplitude"]), q(family["input_bound"])
        scale_interval = _constants(bits)[1]
        scale, scale_radius = _midpoint(scale_interval)
        for i in range(modes):
            k = 2*i+1
            pole_interval = exp_negative(k*k*dt, bits)
            pole, pole_radius = _midpoint(pole_interval)
            rho = pole_interval[1]
            if rho >= 1:
                raise ValueError("Precision cannot separate the slowest heat pole from one; increase bits or dt")
            # The sampled source has exact DC gain 1/k^2 for each mode.
            gain = (1-pole)/(k*k)
            weight = scale/(k*k)
            weight_radius = scale_radius/(k*k)
            state_size = initial+forcing/Fraction(k*k)
            # All-time difference of two scalar filters with tied DC gains.
            state_error = (initial*pole_radius/(1-rho)
                           + forcing*pole_radius/Fraction(k*k)*
                           (1/(1-rho)+(1-pole)/(1-rho)**2))
            coefficient_error = _upper(coefficient_error+weight_radius*state_size+abs(weight)*state_error, bits)
            entries.append([i, i, str(pole)])
            b.append([str(gain)])
            c.append(str(weight))
            center.append("0")
            radius.append(str(initial))
            modal.append({"index": k, "pole_interval": [str(v) for v in pole_interval],
                          "weight_interval": [str(v/Fraction(k*k)) for v in scale_interval]})
        # Convexity over cells centered on the omitted odd integers gives
        # sum_{j>=m} (2j+1)^-p <= (2m)^(1-p)/(2(p-1)).
        tail_initial = scale_interval[1]*initial/Fraction(4*modes)
        tail_forcing = scale_interval[1]*forcing/Fraction(48*modes**3)
        budgets = {"omitted_initial": str(_upper(tail_initial, bits)), "omitted_forcing": str(_upper(tail_forcing, bits)),
                   "retained_coefficient_and_sampling": str(coefficient_error)}
        domain = {"initial_center": center, "initial_radius": radius, "input_bound": str(forcing)}
        semantics = ("Spatial mean of theta_t=theta_xx+u(t) on (0,pi), zero boundary temperatures; "
                     "constant initial field a with |a|<=initial_amplitude, |u(t)|<=input_bound. "
                     "Midpoint model assumes zero-order-held inputs of duration dt. Tail bound itself covers measurable bounded inputs continuously.")
        obstruction = "The unrestricted heat input/output map has infinitely many visible odd poles; restricted-domain truncation errors are recorded separately."
    total = sum(q(value) for value in budgets.values())
    model = {"n": modes, "entries": entries, "B": b, "C": [c]}
    return family, model, domain, modal, budgets, total, semantics, obstruction


def compile_spectral(family, modes, tolerance, bits=80):
    """Create a rational model whose source error is bounded for all future samples."""
    tolerance = _nonnegative(tolerance)
    family, model, domain, modal, budgets, total, semantics, obstruction = _build(family, modes, bits)
    return {"schema": "dynacomp.spectral/1", "family": family, "modes": modes, "bits": bits,
            "model": model, "domain": domain, "modal_enclosures": modal, "error_budget": budgets,
            "all_time_upper": str(total), "tolerance": str(tolerance),
            "status": "tolerance_met" if total <= tolerance else "not_certified",
            "source_semantics": semantics, "exact_realization_obstruction": obstruction,
            "source_reference": PAPER if family["kind"] == "ramanujan_log" else "Analytic spectral formula derived in docs/SPECTRAL.md",
            "claim": "Analytic infinite source -> exact rational midpoint model, all future samples; runtime floating roundoff and physical calibration excluded."}


def replay_spectral(packet):
    """Reconstruct every bound and retained coefficient from the source formula."""
    if packet.get("schema") != "dynacomp.spectral/1":
        raise ValueError("Unknown spectral schema")
    fresh = compile_spectral(packet["family"], packet["modes"], packet["tolerance"], packet["bits"])
    if fresh != packet:
        raise ValueError("Spectral source, enclosure, domain, budget or status differs")
    return True


def search_spectral(family, tolerance, max_modes=256, bits=80):
    """Geometric budgeted search. Success is certified; no rank optimality claim."""
    _integer(max_modes, "Maximum modes")
    modes, history = 1, []
    while True:
        packet = compile_spectral(family, modes, tolerance, bits)
        history.append({"modes": modes, "all_time_upper": packet["all_time_upper"], "status": packet["status"]})
        if packet["status"] == "tolerance_met" or modes == max_modes:
            return packet, history
        modes = min(max_modes, modes*2)


def spectral_error_at(packet, step):
    """A refined source-only upper bound at a specified sample/index."""
    replay_spectral(packet)
    _integer(step, "Step", low=0, high=1000000)
    budget = packet["error_budget"]
    family = packet["family"]
    if family["kind"] == "heat_rod":
        rho = exp_negative((2*packet["modes"]+1)**2*q(family["dt"]), packet["bits"])[1]
    else:
        rho = Fraction(1, (packet["modes"]+1)**2)
    return (q(budget["omitted_initial"])*rho**step+q(budget["omitted_forcing"])
            + q(budget["retained_coefficient_and_sampling"]))


def compress_spectral_modes(packet, keep, tolerance):
    """Keep selected zero-based modal coordinates with an all-time guarantee.

    This uses scalar DC gains, not an arbitrary projected-model stability
    inference. Distinct rational poles mean no exact generic closure is lost
    silently: every deleted contribution is charged to the output budget.
    """
    replay_spectral(packet)
    keep = list(keep)
    if not keep or any(type(i) is not int or not 0 <= i < packet["modes"] for i in keep) or len(set(keep)) != len(keep):
        raise ValueError("Distinct nonempty zero-based modal indices required")
    tolerance = _nonnegative(tolerance)
    model, domain = packet["model"], packet["domain"]
    poles = [q(entry[2]) for entry in model["entries"]]
    c = [q(v) for v in model["C"][0]]
    removed = Fraction(0)
    selected = set(keep)
    for i in range(packet["modes"]):
        if i in selected:
            continue
        size = abs(q(domain["initial_center"][i]))+q(domain["initial_radius"][i])
        if q(domain["input_bound"]):
            size += abs(q(model["B"][i][0]))*q(domain["input_bound"])/(1-poles[i])
        removed = _upper(removed+abs(c[i])*size, packet["bits"])
    e = [[str(int(i == j)) for j in range(packet["modes"])] for i in keep]
    g = [[str(poles[i] if i == j else 0) for j in keep] for i in keep]
    h = [model["B"][i] for i in keep]
    d = [[model["C"][0][i] for i in keep]]
    total = q(packet["all_time_upper"])+removed
    return {"schema": "dynacomp.spectral-subset/1", "source": packet, "keep": keep,
            "E": e, "G": g, "H": h, "D": d, "rank": len(keep),
            "removed_modes_upper": str(removed), "combined_all_time_upper": str(total),
            "tolerance": str(tolerance), "status": "tolerance_met" if total <= tolerance else "not_certified",
            "claim": "Infinite source -> selected rational modes at every future sample; analytic source and deleted-mode budgets added; no arbitrary-projection or rank-optimality claim."}


def replay_spectral_subset(packet):
    if packet.get("schema") != "dynacomp.spectral-subset/1":
        raise ValueError("Unknown spectral subset schema")
    if compress_spectral_modes(packet["source"], packet["keep"], packet["tolerance"]) != packet:
        raise ValueError("Modal selection, all-time budget or acceptance differs")
    return True


def compose_spectral(packet, E, G, H, D, horizon, tolerance, **options):
    """Combine the infinite-source budget with any finite reduced model.

    Finite acceptance uses an enclosing initial box (heat's correlated initial
    field is deliberately enlarged), independent bounded ZOH input samples.
    A failed finite lower witness does not by itself refute the infinite source.
    """
    replay_spectral(packet)
    _integer(horizon, "Horizon", low=0, high=10000)
    tolerance = _nonnegative(tolerance)
    source_error = q(packet["all_time_upper"])
    domain = packet["domain"]
    finite = validated_box(SparseLinearModel.from_spec(packet["model"]), E, G, H, D,
                           domain["initial_center"], domain["initial_radius"],
                           [[domain["input_bound"]]]*horizon,
                           max(Fraction(0), tolerance-source_error), **options)
    combined = source_error+q(finite["certificate"]["worst_upper"])
    return {"schema": "dynacomp.spectral-composition/1", "source": packet, "finite": finite,
            "tolerance": str(tolerance), "combined_upper": str(combined), "horizon": horizon,
            "status": "tolerance_met" if combined <= tolerance and finite["status"] == "tolerance_met" else "not_certified",
            "claim": "Infinite analytic source -> supplied reduction at samples 0..horizon; source-tail and finite-reduction errors added; no all-time claim for second-stage reduction."}


def replay_spectral_composition(packet):
    if packet.get("schema") != "dynacomp.spectral-composition/1":
        raise ValueError("Unknown spectral composition schema")
    replay_spectral(packet["source"])
    replay_validated(packet["finite"])
    finite = packet["finite"]
    fresh = compose_spectral(packet["source"], *(finite[k] for k in ("E", "G", "H", "D")),
                             packet["horizon"], packet["tolerance"],
                             precisions=(finite["certificate"]["bits"],),
                             require_contraction=finite["require_contraction"])
    expected = fresh["finite"]
    # Precision search history is diagnostic; arithmetic/status replayed above.
    fresh["finite"] = finite
    if fresh != packet:
        raise ValueError("Composed domain, source, bound or status differs")
    # Explicitly tie the independently checked finite packet to this source.
    for key in ("model", "initial_center", "initial_radius", "input_radii", "input_centers", "tolerance", "certificate", "status", "rank"):
        if expected[key] != finite[key]:
            raise ValueError("Finite certificate is not for the declared infinite-source domain")
    return True


def simulate_spectral(packet, inputs, initial_amplitude=0, allow_unaccepted=False):
    """Exact rational execution with source-domain and horizon guards.

    initial_amplitude selects the common initial rod field, not arbitrary
    independent modal values. Moment sources have their prescribed initial
    state and require initial_amplitude=0 and every supplied input zero.
    """
    if type(allow_unaccepted) is not bool:
        raise ValueError("Boolean acceptance override required")
    readers = {"dynacomp.spectral/1": replay_spectral,
               "dynacomp.spectral-subset/1": replay_spectral_subset,
               "dynacomp.spectral-composition/1": replay_spectral_composition}
    schema = packet.get("schema")
    if schema not in readers:
        raise ValueError("Unknown spectral execution schema")
    readers[schema](packet)
    if packet["status"] != "tolerance_met" and not allow_unaccepted:
        raise ValueError("Declared tolerance not certified")
    source = packet if schema == "dynacomp.spectral/1" else packet["source"]
    inputs = tuple(q(u) for u in inputs)
    if schema == "dynacomp.spectral-composition/1" and len(inputs) > packet["horizon"]:
        raise ValueError("Execution exceeds certified second-stage horizon")
    bound = q(source["domain"]["input_bound"])
    if any(abs(u) > bound for u in inputs):
        raise ValueError("Input outside infinite-source domain")
    initial = q(initial_amplitude)
    if source["family"]["kind"] == "heat_rod":
        if abs(initial) > q(source["family"]["initial_amplitude"]):
            raise ValueError("Initial field outside infinite-source domain")
        state = (initial,)*source["modes"]
    else:
        if initial:
            raise ValueError("Moment source has a fixed prescribed initial state")
        state = tuple(q(v) for v in source["domain"]["initial_center"])
    if schema == "dynacomp.spectral/1":
        poles = tuple(q(row[2]) for row in source["model"]["entries"])
        b, c = matrix(source["model"]["B"]), matrix(source["model"]["C"])
        def advance(x, u):
            return tuple(a*v+row[0]*u for a, v, row in zip(poles, x, b))
    else:
        reduced = packet["finite"] if schema == "dynacomp.spectral-composition/1" else packet
        state = mv(matrix(reduced["E"]), state)
        g, h, c = (matrix(reduced[key]) for key in ("G", "H", "D"))
        def advance(x, u):
            return addv(mv(g, x), mv(h, (u,)))
    outputs = [mv(c, state)]
    for u in inputs:
        state = advance(state, u)
        outputs.append(mv(c, state))
    return tuple(outputs)
