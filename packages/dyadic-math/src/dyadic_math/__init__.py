from dyadic_math._version import __version__

# 2a. Dynamics
from dyadic_math.basin import BasinExplorer, precision_sweep, LayerGhostDiagnosticV2, GhostHunt
from dyadic_math.separation import (
    newton_trajectory, predicted_separation, verify_separation,
    step_count_profile, ultrametric_ball_tree,
)
from dyadic_math.fourier import (
    step_count_fn, analytic_step_count, dyadic_coefficients,
    fourier_summary, ultrametric_uncertainty,
)
from dyadic_math.mersenne import (
    mersenne_coordinates, mersenne_cliff_table, cliff_constant,
    mersenne_cliff_theorem, cliff_constant_unified, dlog_with_lut,
    prove_cliff_constant, prove_c_formula, proof_connection,
)
from dyadic_math.isometry import (
    verify_isometry, isometry_pair_test, isometry_summary,
)

# 2b. General p-adic
from dyadic_math.padic_roots import (
    lift_root, newton_step, halley_step, newton2_step, newton3_step,
)

# 2c. Weight stability analysis
from dyadic_math.thermodynamics import SeedThermodynamics

# 2d. Extensions
from dyadic_math.crt import CRTDualNumber, CRTDualProcessor, combined_stability
from dyadic_math.nonabelian import NonAbelianCRTDual, phase_alignment_experiment
from dyadic_math.iwasawa import (
    congruence_depth, ldu_decompose, matrix_coordinates,
    matrix_commutator, MatrixCoordinates,
)

__all__ = [
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
