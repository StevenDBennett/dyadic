# Changelog

All notable changes to dyadic are documented here.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added

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

### Changed

- Extracted from the dual-view project (three predecessor codebases)
- Deleted dyadic-ml, dyadic_math.operators, dyadic_math.gauge, dyadic_math.iwasawa_algebra
- Replaced degenerate eps_crit metric with phase_alignment_experiment
- Replaced binary ghost_ratio with graded v₂(e) thermodynamic measure

### Fixed

- α=1 sector ghost attractor: Newton iteration now solves 5^e ≡ -a for targets ≡ 3 mod 4
- NewtonProjector modulus mismatch: modinv computed modulo k-2 (exponent domain) not k
- Average operator O(N²) performance: replaced closure chain with direct summation
- Closure late-binding bug in operator composition
- Test 11 tautology: proper separation of training and test data
