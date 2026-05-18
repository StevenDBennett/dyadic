# Changelog

All notable changes to dyadic are documented here.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added

- Public `dlog_newton_step` function extracted from private `_dlog_newton_step` for reuse in `mersenne.dlog_with_lut`
- `_step` module: `newton_step_core` extracted from `basin.py` to break hub dependency (separation.py, fourier.py now import from `_step`)
- Two-package library structure: dyadic-core (stdlib) and dyadic-math (numpy)
- 2-adic arithmetic: DualNumber decomposition, modular inverse via Newton lifting, discrete logarithms with LUT bootstrap, p-adic exp/log
- Mahler calculus: binomial basis with Dirac/Volterra operators
- Exponent space: additive coordinate chart on Z/2^(k-2)
- Newton dynamics: BasinExplorer, ghost attractor detection, precision sweeps
- Trajectory separation theorem: exact zero-variance verification
- Fourier analysis of step-count functions: analytic vs numeric DFT comparison
- Mersenne Ghost Theorem: cliff constant proofs, bootstrap optimality analysis
- CRT extensions: Z/(2^k · p)Z dual-number system with primitive root dlog
- GL(2) holonomy: non-abelian CRT-dual system with phase alignment experiments
- Iwasawa filtrations: congruence depth, LDU decomposition, commutator depth theorem
- General p-adic root finding: Newton, Halley, composed methods for p ≠ 2, 3
- Seed thermodynamics: graded v₂(e) stability diagnostics, precision-sweep ghost cliff detection
- Weight stability analysis: z-score comparison against random baseline
- Core module split: `dyadic_core.core` → 6 focused submodules (`_exceptions`, `_util`, `_modular`, `_series`, `_dlog`, `_dual`) with backward-compat re-export shim
- Shared `_series_accumulate` helper extracted from duplicated `padic_exp`/`padic_log` loops
- Shared `newton_step_core` helper extracted from basin/fourier/separation duplicate loops
- `__all__` defined in all 6 core submodules
- 11 new tests: `from_coords` edge cases, `pow(0)`, `dlog_bootstrap` direct coverage
- `docs/api.md` dyadic-core section rewritten to match 6-submodule layout
- AGENTS.md: code quality standards section, ruff checks, `__all__` requirement, docs sync verification

### Changed

- `LayerGhostDiagnosticV2`: removed dead `g` and `max_iter` parameters (class only uses `two_adic_dlog`)
- `GhostHunt`: removed dead `max_iter` parameter
- `SeedThermodynamics`: raises `ValueError` if `g != 5` (dlog infrastructure is base-5 only)
- `nonabelian.invariants`: `alpha_det`/`e_det` set to `None` for even determinants (was silently returning coordinates of forced-odd value)
- `ExponentSpace.integrate`: warns instead of raising `ValueError` for large k
- `docs/api.md`: updated `dlog_newton_step`, `LayerGhostDiagnosticV2`, `GhostHunt`, `combined_stability`, `ExponentSpace.integrate` documentation
- `AGENTS.md`: updated test counts (115 core + 108 math), mypy fallback note, constructor validation standard, docs verification procedure
- Extracted from the dual-view project (three predecessor codebases)
- Deleted dyadic-ml, dyadic_math.operators, dyadic_math.gauge, dyadic_math.iwasawa_algebra
- Replaced degenerate eps_crit metric with phase_alignment_experiment
- Replaced binary ghost_ratio with graded v₂(e) thermodynamic measure
- `exp2_neg4(k)` delegates to `g0(k)` instead of duplicate implementation
- `BasinExplorer.newton_step` and `_trajectory` use `self.g` not hardcoded `5`
- `dlog_residual_tracking` bootstrap unified to `k // 2 + 2` (same as `dlog_newton`)
- `mersenne.cliff_constant` for `g != 5` uses `padic_log(g, k)` instead of standalone loop
- `fourier.step_count_fn` and `separation.newton_trajectory` delegate to `newton_step_core`
- Renamed `prove_cliff_constant` → `verify_cliff_constant`, `prove_c_formula` → `verify_c_formula`, `proof_connection` → `verify_connection`
- `isometry.py`: removed T6-b docstring (no implementation exists)
- `iwasawa.py`: `verify_commutator_depth` tests non-commuting shear matrices

### Fixed

- `nonabelian.invariants`: `alpha_det` and `e_det` now correct for odd determinants (was computing on `det_mod2k | 1`)
- `phase_alignment_experiment`, `trace_alpha_independence`, `trace_exponent_independence`: skip trials with even holonomy determinant
- `mersenne.dlog_with_lut`: deduplicated Newton step code (uses `dlog_newton_step` now)
- α=1 sector ghost attractor: Newton iteration now solves 5^e ≡ -a for targets ≡ 3 mod 4
- NewtonProjector modulus mismatch: modinv computed modulo k-2 (exponent domain) not k
- Average operator O(N²) performance: replaced closure chain with direct summation
- Closure late-binding bug in operator composition
- Test 11 tautology: proper separation of training and test data
- `thermodynamics.weight_coordinates` returns `None` for zero (was `float("inf")`)
- `thermodynamics`: fixed `wbm` typo → `w_mod`
- `_dual.from_coords` rejects `v < 0` (was silently accepted)
- `two_adic_log5(10)` gave 901 instead of 389 (pre-existing bug in duplicate `exp2_neg4`)
- `nonabelian.invariants`: narrowed bare `except Exception` to `except (ValueError, ArithmeticError)`
- `pyproject.toml`: removed triple blank line
