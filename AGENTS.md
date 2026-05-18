# dyadic — Agent Guide

Two packages only:

| Package | Location | Description |
|---------|----------|-------------|
| `dyadic-core` | `packages/dyadic-core/` | 2-adic foundation — DualNumber, modinv, dlog, padic exp/log, Mahler calculus, exponent space. Stdlib only |
| `dyadic-math` | `packages/dyadic-math/` | p-adic mathematics: Newton dynamics, CRT extensions, GL(2) holonomy, Iwasawa filtrations, Fourier analysis, Mersenne theorems, general p-adic root finding (p ≠ 2, 3), thermodynamics. Depends on dyadic-core + numpy |

If a module cannot honestly deliver what its name claims, delete it.

## Language

Use **UK English** throughout (code, docs, commit messages). The user lives in Dublin, Ireland.

## Research Code

Code kept alive but unpublished lives in `research/` — not installed, not tested in CI, not referenced in published docs.

## Documentation

- `docs/api.md` documents every public function/class in both packages.
- `docs/bug_history.md` documents known bugs honestly.
- `CHANGELOG.md` describes only what currently exists.
- Do not create `.md` files unless explicitly asked.
- Before committing, verify `docs/api.md` matches the API: every name listed must resolve via `getattr(package, name)`.
- Internal submodules (prefixed `_`) are not public API and must not appear in `docs/api.md`.

## Pre-Commit Checks

Run all three before any commit:
1. `ruff check packages/` — must pass clean
2. `ruff format --check packages/` — must pass clean
3. `python3 -m pytest --tb=short` — 223 tests (115 core + 108 math), all must pass

Mypy `--strict` is configured but unavailable in this environment. If it becomes available, fix all mypy errors too.

## Commit Style

- Describe what changed and why. UK English. No emoji.
- One commit per logical change.

## Code Quality

- Every public submodule must define `__all__`.
- No bare `except Exception` — narrow to `ValueError`, `ArithmeticError`, etc.
- No silent API misbehaviour: reject unsupported parameters early with `ValueError`. Return `None` for inherently undefined values (e.g. valuation of zero), never sentinels like `float("inf")`.
- Do not duplicate iteration loops — extract shared helpers.
- No unused imports.
- If a parameter would silently produce wrong results for most plausible values, reject it in the constructor (e.g. generators other than base-5 in dlog-hardcoded components).

## What Not To Do

- Do not create `.md` files proactively.
- Do not add comments unless they explain a non-obvious invariant.
