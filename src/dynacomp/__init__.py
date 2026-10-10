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

from .robust import switching_enrichment, replay_switching
from .stability import (schur_stable, hurwitz_stable, quadratic_certificate, replay_quadratic,
                        continuous_certificate, replay_continuous, common_quadratic_certificate,
                        replay_common_quadratic)
from .forever import all_time_bound, all_time_enrichment, replay_all_time, moment_rows
from .zonotope import sharp_zonotope_bound, replay_zonotope, box_generators

from .validated import validated_box,replay_validated,BoxKernelCache,euclidean_contraction
from .stable_proposals import stable_observable_proposal,StableProposalContext

from .spectral import (compile_spectral, replay_spectral, search_spectral,
                       compose_spectral, replay_spectral_composition,
                       spectral_error_at, ramanujan_coefficients)
from .spectral import compress_spectral_modes, replay_spectral_subset, simulate_spectral
from .ramanujan import ramanujan_value, replay_ramanujan_value

from .generating import (RationalGF, from_recurrence, from_state, periodic_gf,
                         binomial_tail, replay_generating)
from .machine import (certify_execution, replay_execution, compose_execution,
                      replay_machine_composition, Binary64Runtime)

from .native import (certify_sequential, replay_sequential, compose_sequential,
                     replay_sequential_composition, PreparedSequential)

__version__ = "0.11.0"
