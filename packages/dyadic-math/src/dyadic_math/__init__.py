from dyadic_math._version import __version__

# 2a. Dynamics
from dyadic_math.basin import BasinExplorer, GhostHunt, LayerGhostDiagnosticV2, precision_sweep

# 2d. Extensions
from dyadic_math.crt import CRTDualNumber, CRTDualProcessor, combined_stability
from dyadic_math.fourier import (
    analytic_step_count,
    dyadic_coefficients,
    fourier_summary,
    step_count_fn,
    ultrametric_uncertainty,
)
from dyadic_math.isometry import (
    isometry_pair_test,
    isometry_summary,
    verify_isometry,
)
from dyadic_math.iwasawa import (
    MatrixCoordinates,
    congruence_depth,
    ldu_decompose,
    matrix_commutator,
    matrix_coordinates,
)
from dyadic_math.mersenne import (
    cliff_constant,
    cliff_constant_unified,
    dlog_with_lut,
    mersenne_cliff_table,
    mersenne_cliff_theorem,
    mersenne_coordinates,
    proof_connection,
    prove_c_formula,
    prove_cliff_constant,
)
from dyadic_math.nonabelian import NonAbelianCRTDual, phase_alignment_experiment

# 2b. General p-adic
from dyadic_math.padic_roots import (
    halley_step,
    lift_root,
    newton2_step,
    newton3_step,
    newton_step,
)
from dyadic_math.separation import (
    newton_trajectory,
    predicted_separation,
    step_count_profile,
    ultrametric_ball_tree,
    verify_separation,
)

# 2c. Weight stability analysis
from dyadic_math.thermodynamics import SeedThermodynamics

__all__ = [
    "__version__",
    "BasinExplorer", "precision_sweep", "LayerGhostDiagnosticV2", "GhostHunt",
    "newton_trajectory", "predicted_separation", "verify_separation",
    "step_count_profile", "ultrametric_ball_tree",
    "step_count_fn", "analytic_step_count", "dyadic_coefficients",
    "fourier_summary", "ultrametric_uncertainty",
    "mersenne_coordinates", "mersenne_cliff_table", "cliff_constant",
    "mersenne_cliff_theorem", "cliff_constant_unified", "dlog_with_lut",
    "prove_cliff_constant", "prove_c_formula", "proof_connection",
    "verify_isometry", "isometry_pair_test", "isometry_summary",
    "lift_root", "newton_step", "halley_step", "newton2_step", "newton3_step",
    "SeedThermodynamics",
    "CRTDualNumber", "CRTDualProcessor", "combined_stability",
    "NonAbelianCRTDual", "phase_alignment_experiment",
    "congruence_depth", "ldu_decompose", "matrix_coordinates",
    "matrix_commutator", "MatrixCoordinates",
]
