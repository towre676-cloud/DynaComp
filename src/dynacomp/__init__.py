"""DynaComp: compile only the dynamics needed by your target."""
from .linear import LinearModel, compile_linear, replay_linear, closure_witness, distinguishing_horizon
from .finite import FiniteModel, compile_finite, replay_finite
from .polynomial import (compile_polynomial, replay_polynomial,
                         compile_controlled_polynomial, replay_controlled_polynomial)
from .bounds import finite_horizon_bound, box_horizon_bound, usable_horizon
from .adaptive import adaptive_coordinates, optimal_coordinates, replay_adaptive
from .families import compile_family, replay_family, closure_conditions
from .memory import compile_output_memory, replay_memory
from .realization import minimal_realization, replay_realization
from .switched import (compile_switched, replay_switched, WeightedAutomaton, minimize_automaton,
                       replay_automaton, check_equivalence, replay_equivalence)
from .parametric import parametric_analysis, replay_parametric, classify
from .enrich import greedy_enrichment, replay_enrichment

from .network import Network, compile_network, replay_network, run_scenarios

from .uniform import SparseLinearModel, sharp_box_bound, uniform_enrichment, replay_uniform, output_box_range
from .sparse_family import compile_sparse_family, replay_sparse_family

__version__ = "0.5.0"
