# dyadic — Agent Guide

Two packages only:

| Package | Location | Description |
|---------|----------|-------------|
| `dyadic-core` | `packages/dyadic-core/` | 2-adic foundation — DualNumber, modinv, dlog, padic exp/log, Mahler calculus, exponent space. Stdlib only |
| `dyadic-math` | `packages/dyadic-math/` | p-adic mathematics: Newton dynamics, CRT extensions, GL(2) holonomy, Iwasawa filtrations, Fourier analysis, Mersenne theorems, general p-adic root finding (p ≠ 2, 3), thermodynamics. Depends on dyadic-core + numpy |

If a module's name is actively misleading about its core functionality, rename or delete it.

## Language

Use **UK English** throughout (code, docs, commit messages). The user lives in Dublin, Ireland.

## Research Code

Code kept alive but unpublished lives in `research/` — not installed, not tested in CI, not referenced in published docs.

## Documentation

- `docs/api.md` documents every public function/class in both packages.
- `docs/bug_history.md` documents known bugs honestly.
- `CHANGELOG.md` describes only what currently exists.
- Do not create `.md` files unless explicitly asked.
- Before committing, verify `docs/api.md` matches the API: every name in `__all__` must appear in `docs/api.md`, and every documented public name must resolve via `getattr(package, name)`. Run a quick script to check — don't rely on memory.
- Internal submodules (prefixed `_`) may appear in `docs/api.md` as implementation notes but must be clearly marked as internal. Public API entries must be accessible at the package level via `__all__`.

## Pre-Commit Checks

Run all three before any commit:
1. `ruff check packages/` — must pass clean
2. `ruff format --check packages/` — must pass clean
3. `python3 -m pytest --tb=short` — all tests must pass

Mypy `--strict` is configured but unavailable in this environment. If it becomes available, fix all mypy errors too.

## Commit Style

- Describe what changed and why. UK English. No emoji.
- One commit per logical change.

## Code Quality

- The package `__init__.py` must define `__all__`. Submodules should define `__all__` only if they are a realistic target for `from x import *`.
- No bare `except Exception` — narrow to `ValueError`, `ArithmeticError`, etc.
- No silent API misbehaviour: reject unsupported parameters early with `ValueError`. Return `None` for inherently undefined values (e.g. valuation of zero), never sentinels like `float("inf")`.
- Do not duplicate iteration loops — extract shared helpers.
- No unused imports.
- Tests must be capable of failing. Verifying a known theorem against an implementation is fine (the implementation can be wrong even if the theorem isn't). Do not write tests that pass by construction — e.g. fitting distribution parameters and testing against the same data.
- Docstrings must describe what the code does, not metaphorically adorn it. If a class is called `Thermodynamics` it should genuinely involve thermodynamic concepts; if it computes descriptive statistics, name it accordingly.
- Delete dead code (diagnostic stubs, duplicated self-test functions, research tangents) when found. Do not leave it commented out.
- Public constructors and methods whose cost scales as O(2^k) must document the practical limit and warn (or reject with `ValueError` for truly catastrophic scaling) when invoked near or past it.
- If a parameter would silently produce wrong results for most plausible values, reject it in the constructor (e.g. generators other than base-5 in dlog-hardcoded components).

## What Not To Do

- Do not create `.md` files proactively.
- Do not add comments unless they explain a non-obvious invariant.
