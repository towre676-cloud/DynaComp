"""DynaComp: compile only the dynamics needed by your target."""
from .linear import LinearModel, compile_linear, replay_linear, closure_witness, distinguishing_horizon
from .finite import FiniteModel, compile_finite, replay_finite
from .polynomial import compile_polynomial, replay_polynomial
from .bounds import finite_horizon_bound, usable_horizon

__version__ = "0.1.0"
