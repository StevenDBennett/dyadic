# dyadic — Agent Guide

## Project Identity

2-adic arithmetic library with Newton dynamics analysis. Nothing more.

- **Not** spectral geometry, gauge theory, Iwasawa algebras, ML, quantum compilation, or butterfly compilers.
- If a module cannot honestly deliver what its name claims, delete it — do not keep broken code.

## Package Structure

Two packages only:

| Package | Location | Description |
|---------|----------|-------------|
| `dyadic-core` | `packages/dyadic-core/` | 2-adic arithmetic foundation — stdlib only |
| `dyadic-math` | `packages/dyadic-math/` | Newton dynamics, CRT, Fourier, Mersenne — depends on dyadic-core + numpy |

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

## Commit Style

- Describe what changed and why.
- Use UK English.
- No emoji.
- One commit per logical change.

## What Not To Do

- Do not add spectral geometry, gauge theory, or Iwasawa algebra claims anywhere.
- Do not mention `dyadic-ml`, `butterfly`, `ghost_regularization` as features — they were deleted.
- Do not install or reference `torch` or `scipy` as dependencies.
- Do not create `.md` files proactively unless asked.
- Do not add comments to code unless they explain a non-obvious invariant.
