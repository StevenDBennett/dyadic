> This document is part of the dyadic (2-package library) project, ported from dual-view.

# Research Status

*Research session completed 2026-05-11. Ported from dual-view.*

---

## Completed

### Theorems Proven

| Theorem | Result | Proof |
|---------|--------|-------|
| $Φ_n^*(0)$ | $2^{\mu_3(n)/6}$ | $A_d(0) = 2^{(3^{d-1}-1)/2}$, Möbius inversion |
| $Φ_n^*(1)$ | $3^{\mu_3(n)/2}$ | $B_d(1) = 3^{(3^d-1)/2}$, L'Hôpital |
| $\prod\mu$ | $6^{\mu_3(n)/6}$ | Ratio of $Φ_n^*(1)/Φ_n^*(0)$ |
| $\mu_u = \mu_x$ | Period-4 identity | Cycle equation $\prod(2u_i+1) = 81\prod u_i$ |

### Computations

- **Period-5 dynatomic**: Degree 80 in $u$, irreducible, discriminant 5,605 digits
- **Period-6 dynatomic**: Degree 232 in $u$, coefficients computed (233 entries)
- **Clean primes**: Confirmed {7, 103, 181} below 100,000
- **Galois group (period-4)**: Subgroup of $S_6 \wr V_4$

| Period | deg(u) | $\Phi(0)$ | $\Phi(1)$ | $\prod\mu$ | $\sum\mu$ | Cycles |
|--------|--------|-----------|-----------|------------|-----------|--------|
| 2 | 2 | $2^1$ | $3^3$ | $6^1$ | 6 | 1 |
| 4 | 24 | $2^{12}$ | $3^{36}$ | $6^{12}$ | 90 | 6 |
| 5 | 80 | $2^{40}$ | $3^{120}$ | $6^{40}$ | 486 | 16 |
| 6 | 232 | $2^{116}$ | $3^{348}$ | $6^{116}$ | — | — |

## In Progress

### Galois Group Identification
- **Evidence**: Subgroup of $S_6 \wr V_4$, order 2,949,120
- **Next**: Check parity of block permutations; feed constraints into GAP/Magma for exact identification

### Period-6 Dynatomic
- **Status**: Coefficients computed (degree 232 in $u$)
- **Strategy**: Evaluation-interpolation to avoid coefficient explosion
- **Next**: Multiplier computation from coefficients

## Open Problems

### Critical
1. **Identify exact period-4 Galois group** — feed constraints into GAP/Magma
2. **Prove clean primes = {7, 103, 181}** — use obstruction accumulation
3. **Prove $Φ_n^*(1) = 3^{\mu_3(n)/2}$** — full theoretical proof (currently has complete proof — see `theorems.md`)

### Known Caveats
1. **Fast-check false positives**: $p = 313$ passes periods 2–5 but has a period-27 ghost cycle. Full clean-prime verification requires checking **all** cycles and multipliers, not just root existence modulo $p$.
2. **$p = 181$ cube-resolution**: The period-4 polynomial has 4 linear factors mod 181, but none are cubes in $\mathbb{F}_{181}$ — without the cube-root check, this would appear to have period-4 points. The cube-root test is essential.
3. **3-adic quadratic convergence**: $N(1+\varepsilon) \approx 1 + \varepsilon^2$ implies $v_3(N(1+\varepsilon)-1) = 2 \cdot v_3(\varepsilon)$ — the Newton map is quadratically convergent in the 3-adic metric, consistent with its superattracting fixed points at cube roots of unity.

### High
4. **Compute period-6 multipliers** — numerical computation from coefficients
5. **Find closed form for individual multipliers** — not just product

### Medium
6. **Relate $K_2$ to ray class fields of $\mathbb{Q}(\zeta_3)$**
7. **Generalize to degree $d$ maps**
