> This document is part of the dyadic (3-package monorepo) project, ported from dual-view.

# Computation Data

## Degree Formula

$$\mu_3(n) = \sum_{d|n} \mu(n/d) \cdot 3^d$$

| Period | $\mu_3(n)$ | deg(x) | deg(u) |
|--------|------------|--------|--------|
| 2 | 6 | 6 | 2 |
| 3 | 24 | 24 | 8 |
| 4 | 72 | 72 | 24 |
| 5 | 240 | 240 | 80 |
| 6 | 696 | 696 | 232 |

## Special Values

| Period | $\Phi_n^*(0)$ | $\Phi_n^*(1)$ | $\prod\mu$ |
|--------|---------------|---------------|------------|
| 2 | $2^1 = 2$ | $3^3 = 27$ | $6^1 = 6$ |
| 4 | $2^{12} = 4,\!096$ | $3^{36}$ | $6^{12} = 2,\!176,\!782,\!336$ |
| 5 | $2^{40}$ | $3^{120}$ | $6^{40}$ |
| 6 | $2^{116}$ | $3^{348}$ | $6^{116}$ |

## Period 2

- **Polynomial**: $20u^2 + 5u + 2$ in $u = x^3$
- **Galois group**: $D_6$ (dihedral, order 12), field $K_2 = \mathbb{Q}(\zeta_3, \sqrt{5}, u_1^{1/3})$ where $u_1 = (-5 + 3\sqrt{-15})/40$, $[K_2:\mathbb{Q}] = 12$
- **Multiplier**: $\mu = 6$ (universal), $\sum\mu = 6$

## Period 3

- **Polynomial**: $(19u^2 + 7u + 1) \cdot P_{18}(u)$ in $u = x^3$
- **First factor Galois group**: $C_3 \wr C_2$ (wreath product, order 18), field $K_3 = \mathbb{Q}(\zeta_3, v_1^{1/3}, 19^{1/3})$ where $v_1 = (-7 + 3\sqrt{-3})/38$, $[K_3:\mathbb{Q}] = 18$
- **Multiplier**: $\mu = \pm 24\sqrt{-3}$, $|\mu| = 24\sqrt{3} \approx 41.57$

## Period 4

- **Polynomial**: Degree 24 in $u = x^3$
- **Irreducible**: Yes, over $\mathbb{Q}$
- **Discriminant**: 506 digits, **not a square**
- **Leading coefficient**: $2^{34} \cdot 13 = 223,\!338,\!299,\!392$
- **Constant term**: $2^{12} = 4,\!096$
- **Galois group**: Subgroup of $S_6 \wr V_4$ (wreath product, order 2,949,120), $[K_4:\mathbb{Q}] = 24$
- **Linear disjointness**: The period-4 polynomial is irreducible over $\mathbb{Q}(\zeta_3)$, $\mathbb{Q}(\sqrt{5})$, $\mathbb{Q}(\sqrt{-15})$, $\mathbb{Q}(\zeta_3, \sqrt{5})$ (the period-2 base field), and $\mathbb{Q}(\zeta_3, 19^{1/3})$ (the period-3 base field) — implying $K_4$ is linearly disjoint from $K_2$ and $K_3$ over $\mathbb{Q}(\zeta_3)$

### Multipliers

Minimal polynomial:
$$\mu^6 - 90\mu^5 + 1134\mu^4 - 431568\mu^3 + 18475776\mu^2 - 287214336\mu + 2176782336 = 0$$

| Pair | $\mu$ | $|\mu|$ |
|------|-------|---------|
| 1,2 | $-28.25 \mp 59.94i$ | 66.27 |
| 3,4 | $9.48 \pm 10.83i$ | 14.39 |
| 5 | $22.87$ | 22.87 |
| 6 | $104.67$ | 104.67 |

$\sum\mu = 90$, $\prod\mu = 6^{12}$

### Galois evidence
- No 24-cycles found in 2,261 primes tested (Chebotarev)
- Fixed points always multiples of 4
- 24 roots form 6 blocks of 4, matching $S_6 \wr V_4$

## Period 5

- **Polynomial**: Degree 80 in $u = x^3$, irreducible over $\mathbb{Q}$
- **Discriminant**: 5,605 digits, **not a square**
- **Leading coefficient**: $2^{116} \cdot 211$
- **Constant term**: $2^{40}$

### Multipliers (16 cycles)

| Pair | $\mu$ | $|\mu|$ |
|------|-------|---------|
| 1,2 | $9.11 \pm 16.17i$ | 18.56 |
| 3,4 | $33.94 \pm 11.06i$ | 35.70 |
| 5,6 | $55.87 \pm 33.31i$ | 65.05 |
| 7,8 | $-49.65 \mp 65.21i$ | 81.96 |
| 9,10 | $-94.74 \mp 68.24i$ | 116.76 |
| 11,12 | $118.51 \mp 83.05i$ | 144.71 |
| 13,14 | $188.86 \pm 37.20i$ | 192.49 |
| 15,16 | $-18.89 \mp 317.63i$ | 318.19 |

$\sum\mu = 486 = 2 \cdot 3^5$, $\prod\mu = 6^{40}$

All $|\mu| > 1$ — all period-5 points are ghost cycles.

## Period 6

- **Polynomial**: Degree 232 in $u = x^3$ (233 coefficients, computed)
- **Predicted**: $\Phi(0) = 2^{116}$, $\Phi(1) = 3^{348}$, $\prod\mu = 6^{116}$
- **Multipliers**: Not yet computed (coefficients loaded from `data/period6_coefficients_full.txt`)
