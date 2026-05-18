# Bug Fix History and Audit

> This is a dyadic (2-package library) document, ported from the dual-view project.

## Summary

Across the three predecessor projects, **8 bugs** were identified and fixed. The unified dyadic library incorporates all fixes. Two additional bugs (numpy.int64 crash in GhostMap, training history overwrite) affected code that was removed in the refactor and are no longer applicable.

## Critical Bugs

### 1. α=1 Sector Ghost Attractor (BasinExplorer.newton_step)

**Location**: `dyadic_math.basin` / `BasinExplorer`

**Bug**: For targets with `α = 1` (`a ≡ 3 mod 4`), the Newton iteration solved `5^e ≡ a` instead of `5^e ≡ -a`. This caused convergence to a ghost fixed point at `e* = dlog(a+2, k)` for all α=1 targets.

**Fix**: Added check in `__init__`:
```python
if alpha:
    self.a = (-target) & self.mask
else:
    self.a = target
```

**Impact**: Affected v2 and v3 — both lacked this fix. Without it, ALL weights in the α=1 sector (roughly half of all odd weights) would be classified as ghosts, making the GhostMap binary rather than graded.

**Found**: Original v1 codebase contained the fix; v2/v3 regressed when code was rewritten.

### 2. NewtonProjector Modulus Mismatch

**Location**: `dyadic_math.operators` / `NewtonProjector._step`

**Bug**: In v2, the modular inverse was computed modulo `k` (full ring) instead of `k-2` (exponent domain): `modinv_newton(df_unit, self.ctx.k)`. Since `df_unit` is masked to `k-2` bits, passing `k` as the modulus was incorrect.

**Fix**: Changed to `modinv_newton(df_unit, self.ctx.k - 2)`.

**Impact**: Only affected v2. V1 and v3 had the correct modulus.

### 3. OperatorContext.g Attribute Missing

**Location**: `dyadic_math.operators` / `OperatorContext`

**Bug**: The `g` (generator) attribute was computed inside the initialization but not stored as `self.g`, causing `AttributeError` when `NewtonProjector` accessed `ctx.g`.

**Fix**: Added `self._g = g` and `@property def g(self): return self._g`.

**Impact**: v1 original had this bug; fixed in v1 and preserved in v2/v3.

## High-Severity Bugs

### 4. Average Operator O(N²) Performance

**Location**: `dyadic_math.operators` / `OperatorContext.avg`

**Bug**: The average operator was implemented as `(1/N) Σ Sⁱ`, where the `_Operator.__pow__` chained `N` shift operators, each adding a closure layer. For `N = 2^(k-2)`, this created `O(N²)` closures at construction time, which was fatal for `k ≥ 8` (256 closures → 65,536 compositions).

**Fix**: Replaced with a direct O(N) summation lambda:
```python
def avg_action(f, e, ctx=self):
    total = 0
    for t in range(ctx.N):
        total += f(t)
    return total & ctx.mask
```

### 5. Deprecated eps_crit Metric in ramp_break_strength

**Location**: `dyadic_math.nonabelian` / `ramp_break_strength`

**Bug**: The function measured `eps_crit` — the smallest `ε` such that a perturbation `d → d + ε` makes `det` even. For the shear matrix `[[d, 1], [0, 1]]`, the determinant is always `d`, and the smallest `ε` giving an even result is always `1`. The metric is degenerate.

**Fix**: Added `warnings.warn(..., DeprecationWarning)` and replaced with `phase_alignment_experiment()` which measures whether the α-sector of the determinant flips under perturbation.

## Medium-Severity Bugs

### 6. weight_product mod=1 Bug

**Location**: `dyadic_math.gauge` / `weight_product` (v1) — module removed in refactor

**Bug**: For plain integer lists with `mod=1`, the product returned 0 instead of the modular reduction.

**Fix**: Added explicit mod check: `mod = max(mod, 1 << k)`.

### 7. Closures in _Operator.__add__ / __sub__

**Location**: `dyadic_math.operators` / `_Operator.__add__`, `__sub__` — module removed in refactor

**Bug**: Closures captured loop variables by reference (late binding), so all composed operators used the final values.

**Fix**: Added default-argument capture: `def action(f, e, s=self, o=other):`.

### 8. Test 11 Was a Tautology

**Location**: Original test suite (pre-dating dual-view)

**Bug**: "Test 11" tested Newton correction uniformity by comparing the first-step correction distribution to a uniform distribution — but the test used the same data to fit the distribution parameters and test against them, making it a tautology.

**Fix**: Proper separation of training and test data.

## Audit Findings (from v1/v2 Audits)

### Boundary Testing (33 pushes)

The boundary testing campaign found:
- **2 new theorems**: global convergence for α=0, ghost fixed point formula
- **4 bug classes**: all fixed
- **3 fatal flaw analyses**: ghost regularization signal, ramp break metric, popcount compression correlation

### Fatal Flaw: Ghost Regularization Signal

The ghost regularization signal (from the now-removed `GhostMap`) is **binary** after the α=1 fix — every odd weight has convergence ratio 1.0. The genuinely graded stability measure is `v₂(e_true)`, which is the quantity tracked by `SeedThermodynamics`. The original ghost penalty provided no useful gradient for odd weights.

This is why the approach was abandoned in favour of the thermodynamic `v₂(e)` metric. The GhostMap and ghost_penalty modules were removed in the refactor.

### Phase Alternation in GL(2) Holonomy

The α-sector of the holonomy determinant flips under single-bit perturbations at a rate of ~68% (empirical, 30 trials). This phase alignment signal is the replacement for the degenerate `eps_crit` metric.

## Bug Fix Chronology

| # | Bug | Found | Fixed in | Type |
|---|-----|-------|----------|------|
| 1 | α=1 ghost attractor | v1 original | v1 → dyadic | Critical |
| 2 | NewtonProjector modulus | v2 | v1/v3 → dyadic | Critical |
| 3 | ctx.g missing | v1 | v1 → dyadic | Critical |
| 4 | avg O(N²) | v3 addendum | v3 → dyadic | High |
| 5 | eps_crit degenerate | v3 addendum | v3 → dyadic | High |
| 6 | mod=1 bug | v1 | v1 → dyadic | Medium |
| 7 | Closure capture | v1 → v3 | v3 → dyadic | Medium |
| 8 | Test 11 tautology | v3 addendum | v3 audit | Low |

All 8 bugs are fixed in the dyadic library.
