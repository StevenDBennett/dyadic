# Precision Windows in 2-Adic Arithmetic

Unified treatment of the precision-window phenomenon across dyadic, dual-view,
and BigInteger/carryline projects.

## 1. What Is a Precision Window?

2-adic arithmetic with finite word size W exhibits a *precision window*:
intermediate values during ghost-map recovery, basis conversion, or series
truncation require more headroom than the final result. When the intermediate
exceeds `2^W`, truncation corrupts the result irrecoverably.

All precision windows in this ecosystem are instances of a single principle:

> **If an operation divides by `n` after accumulation, the accumulated sum must
> be < `n · 2^W` for exact recovery.**

The three concrete windows below follow from this rule.

---

## 2. Witt Vector Ghost Recovery

Recovery of Witt components from combined ghost values:

```
r_0 = G_0                (mod 2^W)
r_j = (G_j - S_j) / 2^j  (mod 2^W),   where
S_j = Σ_{i<j} 2^i · r_i^{2^{j-i}}
```

The division by `2^j` means `G_j - S_j` must be exactly divisible by `2^j`.
Since `r_i` is stored at word precision `W`, the term `r_i^{2^{j-i}}` occupies
up to `2^{j-i} · (W - i)` bits. After shifting by `2^i`, the sum `S_j` can
require up to `W - j` bits per component. The recovery produces `r_j` with
`W - j` meaningful bits — the **precision window**:

```
r_j < 2^{W-j}
```

When this holds, the division by `2^j` is lossless. When it fails (e.g., a
Witt component is large enough that `G_j - S_j` overflows `gw_t`), the
recovered component is corrupted.

**Reference**: `dyadic.h:1389–1390`, `ghost_recover` at `dyadic.h:1417–1445`.

---

## 3. Taylor Basis Conversion

Basis conversion `Monomial → Taylor → Monomial` uses Stirling numbers:

```
T_j = Σ_{k} s(k, j) · p_k         (monomial → Taylor)
p_i = Σ_{j} S(j, i) · T_j         (Taylor → monomial)
```

The roundtrip requires `T_j · j! < 2^W` for all `j`. The `j!` factor grows
quickly: for `W = 64`, `j!` overflows at `j = 21` (20! < 2^64 < 21!).
When `T_j` is the maximum `FF_j` coefficient, the constraint is:

```
T_j = j! · FF_j < 2^W
```

This is a **factorial precision window**: basis roundtrips are exact only for
polynomials with small enough Taylor coefficients.

**Reference**: `dyadic.h:568–571`, `check_taylor_roundtrip_precision` at
`dyadic.h:1643–1656`.

---

## 4. p-Adic exp/log Term Truncation

The exponential series `exp(x) = Σ x^n/n!` and logarithm
`log(1+y) = Σ (-1)^{n+1} y^n/n` are truncated at `max_terms = 8·sizeof(T)`
(16 for `uint128_t`, the Witt ghost type for `W = uint32_t`).

For `exp(x)` to converge within this truncation, we need:

```
v₂(x^n/n!) = n·v₂(x) - v₂(n!) ≥ 8·sizeof(T)
```

for all remaining terms. The smallest `n` where this holds depends on `v₂(x)`:

| v₂(x) | Terms needed (128-bit) | Status with 16 max |
|-------|----------------------|--------------------|
| ≥ 9   | ≤ 15                 | Full convergence   |
| 8     | ~18                  | Missing 2 terms    |
| 7     | ~21                  | Missing 5 terms    |
| 5     | ~31                  | Missing 15 terms   |
| 3     | ~42                  | Missing 26 terms   |
| 2     | ~63                  | Missing 47 terms   |
| 1     | ~127                 | Missing 111 terms  |

The Witt exp/log roundtrip at 128-bit ghost precision requires:

```
v₂(ghost_j(a)) ≥ 9    for all j
```

For Teichmüller inputs (`a_j = 0` for `j > 0`): `v₂(ghost_j(a)) = 2^j · v₂(a_0)`,
so `v₂(a_0) ≥ 9` suffices. For non-zero higher components, the constraints
tighten: `v₂(a₁) ≥ 8`, `v₂(a₂) ≥ 7`, etc.

**Reference**: `dyadic.h:1765–1777` (`p_adic_exp_impl`), `dyadic.h:1738–1760`
(`p_adic_log_impl`), `dyadic_verify.h:754–768` (proof constraints).

---

## 5. Mersenne Ghost Theorem Window

The Mersenne Ghost Theorem (dual-view) governs the quantization cliff for
`5^k mod 2^n`:

```
5^{2^{n-2}} ≡ 1 - 2^n  (mod 2^{n+1})
```

The **cliff formula** gives the precision boundary:

```
k* = n + 2 + max(0, v₂(n) - 1)
```

The `max(0, v₂(n) - 1)` term is a secondary correction: at `n = 2^m`,
the core identity gains an extra term `2^{2n-4}` which shifts the cliff
when `2n - 4 < n + 2` (i.e. `n < 6`):

```
ε(n) = v₂(n) - 1
```

For n = 4: ε(4) = 1, giving k* = 4 + 2 + 1 = 7 (instead of 6).

**Cliff density**: the expected cliff cost is only 0.25 extra bits,
and Pr[c ≥ 4] ≈ 1.56%.

**Reference**: `dual-view/docs/mersenne_ghost_theorem.md:76–90`.

---

## 6. Bootstrap Optimality

The bootstrap recurrence for computing `5^{2^{n-2}}` step by step:

```
eprec_i = min(2^i · eprec_0, k - 2)
```

Optimal initial precision:

```
eprec_0 = k/2
```

Verified by exhaustive search for all `k ≤ 128`. At `eprec_0 < k/2`, the
bootstrap undershoots and the final result is incorrect; at `eprec_0 > k/2`,
the extra bits are unused.

**Reference**: `dual-view/docs/mersenne_ghost_theorem.md:104–114`.

---

## 7. Summary: Window Table

| Context              | Constraint                           | Critical factor      |
|----------------------|--------------------------------------|----------------------|
| Witt ghost recovery  | `r_j < 2^{W-j}`                     | Component index j    |
| Taylor basis         | `j! · FF_j < 2^W`                   | Degree j             |
| p-adic exp           | `v₂(x) ≥ 8·sizeof(T) / 8`          | Valuation of input   |
| Mersenne cliff       | `k* = n + 2 + max(0, v₂(n)-1)`     | Target bit n         |
| Bootstrap optimality | `eprec_0 = k/2`                     | Available bits k     |

All five are facets of the same underlying constraint: **2-adic arithmetic with
W-bit words has `log₂(W)` bits of headroom for 2-adic operations**, and every
operation consumes a predictable amount of that headroom.

## 8. Mitigations

1. **Widen arithmetic** — Use `gw_t = widen_t<dword_t<W>>` (up to 4× word
   width) for ghost-map accumulation. Already implemented in dyadic.

2. **Dynamic term count** — Compute `max_terms` based on `v₂(x)` instead of
   using a fixed `8·sizeof(T)`. Not yet implemented; would extend Witt exp/log
   coverage to non-Teichmüller inputs.

3. **Increase word size** — Use `W = uint64_t` (128-bit ghosts) instead of
   `W = uint32_t` for more headroom. The `detail::uint128_t` software type
   makes this portable.

4. **Component scaling** — For Witt vectors, constrain higher components:
   `v₂(a_j) ≥ ceil(log₂(N)) - j` ensures ghost-map precision is sufficient
   for Newton recovery.
