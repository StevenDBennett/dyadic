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

## Testing

- Tests live in `packages/dyadic-core/tests/` and `packages/dyadic-math/tests/`.
- Run with: `python3 -m pytest`
- Full suite must pass before any commit: 212 tests (104 core + 108 math) as of May 2026.
- Use `--tb=short` for concise output.

## Workflow

- Create a feature branch for each logical chunk of work (`git checkout -b topic/description`).
- Commit early and often on the branch. One commit per logical change at merge time.
- Squash or rebase before merging to `main` if the branch has many small commits.
- Always run the full test suite (`python3 -m pytest --tb=short`) before merging.

## Commit Style

- Describe what changed and why.
- Use UK English.
- No emoji.
- One commit per logical change.

## What Not To Do

- Do not mention `dyadic-ml`, `butterfly`, `ghost_regularization` as features — they were deleted.
- Do not install or reference `torch` or `scipy` as dependencies.
- Do not create `.md` files proactively unless asked.
- Do not add comments to code unless they explain a non-obvious invariant.
