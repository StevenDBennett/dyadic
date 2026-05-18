# dyadic — Agent Guide

## Project Identity

p-adic arithmetic library with 2-adic Newton dynamics analysis.

- **dyadic-core**: 2-adic arithmetic foundation — DualNumber decomposition, modular inverse via Newton lifting, discrete logarithms, p-adic exp/log, Mahler calculus, exponent space. Stdlib only.
- **dyadic-math**: p-adic mathematics built on top — Newton dynamics (basin analysis, ghost detection), CRT extensions, GL(2) holonomy, Iwasawa-style matrix filtrations, Fourier analysis of step-count functions, Mersenne theorems, general p-adic root finding (Newton, Halley, composed methods for p ≠ 2, 3), weight thermodynamics. Depends on dyadic-core + numpy.

If a module cannot honestly deliver what its name claims, delete it — do not keep broken code.

## Package Structure

Two packages only:

| Package | Location | Description |
|---------|----------|-------------|
| `dyadic-core` | `packages/dyadic-core/` | 2-adic arithmetic foundation — DualNumber, modinv, dlog, padic exp/log, Mahler calculus, exponent space. Stdlib only |
| `dyadic-math` | `packages/dyadic-math/` | p-adic mathematics: Newton dynamics, CRT extensions, GL(2) holonomy, Iwasawa filtrations, Fourier analysis, Mersenne theorems, general p-adic root finding (p ≠ 2, 3), thermodynamics. Depends on dyadic-core + numpy |

`dyadic-ml` was deleted. `dyadic_math.operators`, `dyadic_math.gauge`, `dyadic_math.iwasawa_algebra` were deleted from dyadic-math.

## Language

Use **UK English** throughout (code, docs, comments, commit messages):
- `analyse` not `analyze`
- `colour` not `color` (if applicable)
- `centre` not `center` (if applicable)
- `behaviour` not `behavior` (if applicable)
- `modelled` not `modeled` (if applicable)
- `regularisation` not `regularization` (if applicable)

The user lives in Dublin, Ireland.

## Research Code

Code that is not part of the published packages but is worth keeping lives in `research/`:
- `research/newton_dynamics/` — p-adic Newton dynamics for cube roots (dynatomic polynomials, clean primes)
- `research/unified_butterfly_guidance.md` — butterfly bridge documentation

Research code is not installed, not tested in CI, and not referenced in published docs.

## Documentation Standards

- All docs files describe the current 2-package library only. Never reference dyadic-ml, operators, gauge, iwasawa_algebra, or newton_dynamics (as a subpackage).
- `docs/api.md` documents every public function/class in both packages.
- `docs/bug_history.md` documents known bugs honestly. Entries for deleted code are removed.
- `CHANGELOG.md` describes only what exists.
- Do not create README.md or documentation files unless explicitly asked.
- Before committing, verify `docs/api.md` matches the actual API — every documented top-level name must resolve via `getattr(package, name)`. Run `python3 -c "import dyadic_core, dyadic_math; [hasattr(dyadic_core, n) or hasattr(dyadic_math, n) or print(f'MISSING: {n}') for n in ['name1', ...]]"` or consult the verification script in the agent.
- Internal submodules (prefixed `_`) are not public API and should not appear in `docs/api.md`.

## Testing

- Tests live in `packages/dyadic-core/tests/` and `packages/dyadic-math/tests/`.
- Run with: `python3 -m pytest` or `pytest`.
- Full suite must pass before any commit: 223 tests (115 core + 108 math).
- Use `--tb=short` for concise output.
- Run `ruff check packages/` and `ruff format --check packages/` — both must pass clean before any commit.

## Type Checking

- All source packages are checked with mypy `--strict` (configured in root `pyproject.toml`).
- Test directories and `research/` are excluded from strict checking.
- Run with: `make typecheck` or `python3 -m mypy packages/`
- Test files are not checked in CI (too many untyped test methods).
- Fix all mypy errors before committing source code. If mypy is unavailable in the local environment, at minimum ensure ruff and pytest pass clean.

## Workflow

- Create a feature branch for each logical chunk of work (`git checkout -b topic/description`).
- Commit early and often on the branch. One commit per logical change at merge time.
- Squash or rebase before merging to `main` if the branch has many small commits.
- Before merging, run all three checks: `ruff check packages/`, `ruff format --check packages/`, and `python3 -m pytest --tb=short`. Fix all mypy errors too (or at minimum ensure ruff and pytest pass clean if mypy is unavailable).
- Delete the feature branch after merging.

## Commit Style

- Describe what changed and why.
- Use UK English.
- No emoji.
- One commit per logical change.

## Code Quality

- Every public submodule must define `__all__` listing its exported names.
- Never use bare `except Exception` — narrow to specific exception types (`ValueError`, `ArithmeticError`, etc.).
- No silent API misbehaviour: if a parameter is accepted but ignored, raise `ValueError` or remove it. If a return value is inherently undefined (e.g. valuation of zero), return `None` not a sentinel like `float("inf")`.
- Do not duplicate iteration loops. Extract shared helpers for repeated patterns (series accumulation, Newton step core, etc.).
- No unused imports. MyPy's `--strict` and Ruff will flag these.
- Accepting a parameter that silently produces wrong results for most values is worse than rejecting invalid parameters early. Constructor validation should raise `ValueError` with a clear message when a parameter is unsupported, even if the parameter seems plausibly general (e.g. generator `g` in components that hardcode base-5 dlog).

## What Not To Do

- Do not mention `dyadic-ml`, `butterfly`, `ghost_regularization` as features — they were deleted.
- Do not install or reference `torch` or `scipy` as dependencies.
- Do not create `.md` files proactively unless asked.
- Do not add comments to code unless they explain a non-obvious invariant.
