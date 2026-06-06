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

**Reference**: `include/dyadic/witt.h`, `ghost_recover` function.

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

**Reference**: `include/dyadic/basis.h`, `check_taylor_roundtrip_precision`
(in `include/dyadic/verify.h`).

---

## 4. p-Adic exp/log Term Truncation

The exponential series `exp(x) = Σ x^n/n!` and logarithm
`log(1+y) = Σ (-1)^{n+1} y^n/n` terminate early when term valuation reaches
the bit width of type `T` (128 for `uint128_t`, the Witt ghost type).
A generous budget of `2·bits` terms (256 for 128-bit, 128 for 64-bit) covers
all inputs that can mathematically converge.

Term valuation for exp (`v₂(x^n/n!) = n·v₂(x) - v₂(n!)`) grows linearly when
`v₂(x) ≥ 2`; at `v₂(x) = 1` it never reaches full precision (stalls at
`≤ log₂(n)+1`). The budget scales to the type width at 2×.

| v₂(x) | Terms needed (128-bit) | Budget (2×) |
|-------|----------------------|--------------|
| ≥ 9   | ≤ 15                 | Full         |
| 8     | ~18                  | Full         |
| 5     | ~31                  | Full         |
| 3     | ~42                  | Full         |
| 2     | ~63                  | Full         |
| 1     | ∞ (series stalls)    | Never converges |

For log (`v₂(y^n/n) = n·v₂(y) - v₂(n)`): with `v₂(y) = 1`, needs ~135 terms
at 128-bit (within 256 budget). With `v₂(y) ≥ 2`, needs ~67 terms.

The Witt exp/log roundtrip at 128-bit ghost precision converges for any
`v₂(ghost_j(a)) ≥ 2` (the exp convergence threshold). For `v₂ = 1` inputs,
exp can't converge at any finite precision — this is a mathematical limit,
not a code bug.

For Teichmüller inputs (`a_j = 0` for `j > 0`): `v₂(ghost_j(a)) = 2^j · v₂(a_0)`,
so `v₂(a_0) ≥ 2` suffices for convergence.

**Reference**: `include/dyadic/witt.h` (`p_adic_exp_impl`, `p_adic_log_1plus_impl`),
`include/dyadic/verify.h` (proof constraints).

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

**Reference**: `dual-view/docs/mersenne_ghost_theorem.md`.

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

**Reference**: `dual-view/docs/mersenne_ghost_theorem.md`.

---

## 7. Summary: Window Table

| Context              | Constraint                           | Critical factor      |
|----------------------|--------------------------------------|----------------------|
| Witt ghost recovery  | `r_j < 2^{W-j}`                     | Component index j    |
| Taylor basis         | `j! · FF_j < 2^W`                   | Degree j             |
| p-adic exp           | `v₂(x) ≥ 2, budget = 2·bits`      | Valuation of input   |
| Mersenne cliff       | `k* = n + 2 + max(0, v₂(n)-1)`     | Target bit n         |
| Bootstrap optimality | `eprec_0 = k/2`                     | Available bits k     |

All five are facets of the same underlying constraint: **2-adic arithmetic with
W-bit words has `log₂(W)` bits of headroom for 2-adic operations**, and every
operation consumes a predictable amount of that headroom.

## 8. Mitigations

1. **Widen arithmetic** — Use `gw_t = widen_t<dword_t<W>>` (up to 4× word
   width) for ghost-map accumulation. Already implemented in dyadic.

2. **Increase word size** — Use `W = uint64_t` (128-bit ghosts) instead of
   `W = uint32_t` for more headroom. The `detail::uint128_t` software type
   makes this portable.

3. **Component scaling** — For Witt vectors, constrain higher components:
   `v₂(a_j) ≥ ceil(log₂(N)) - j` ensures ghost-map precision is sufficient
   for Newton recovery.
