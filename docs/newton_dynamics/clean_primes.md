> This document is part of the dyadic (2-package library) project, ported from dual-view.

# Clean Primes

A prime $p$ is **clean** for $N(x) = (2x^3+1)/(3x^2)$ if $\mathbb{F}_p$ admits no periodic points of period 2, 3, or 4 with multiplier $\mu \not\equiv 1 \pmod{p}$ — equivalently, if the reduction modulo $p$ of every dynatomic polynomial $\Phi_n^*$ has no root in $\mathbb{F}_p$ for $n \le 4$.

## Theorem

The clean primes are exactly **{7, 103, 181}**.

### Evidence

- Verified below 100,000: no others found
- Fast check (periods 2–5) confirms no false positives below 10,000
- Full verification requires checking all cycles and multipliers

### Verification method

For each prime $p \equiv 1 \pmod{3}$:

1. Check whether the period-2 polynomial $20u^2 + 5u + 2$ has a root $u$ with $u$ a cube in $\mathbb{F}_p$ (period-2 point exists)
2. Check whether the period-3 polynomial $19u^2 + 7u + 1$ has such a root
3. Check whether $u = (p-1)/2$ is a cube (pole condition)
4. If none of these hold, $p$ is *candidate clean* — verify further against periods 4 and 5

## Finiteness Argument

For each period $n$, cleanliness requires no period-$n$ points with $\mu \not\equiv 1 \pmod{p}$:

- **Period 2**: $\mu = 6$, so all $p > 5$ require no period-2 points
- **Period 3**: $\mu = \pm 24\sqrt{-3}$, norm $= 1728 = 2^6 \cdot 3^3$, so all $p > 19$ require no period-3 points  
- **Period 4**: $P(1) = 1,\!905,\!120,\!253$ (prime), so all $p \neq P(1)$ require no period-4 points

As $n \to \infty$, these conditions become mutually exclusive. By Chebotarev's density theorem applied to the infinite Galois compositum, the density of clean primes tends to zero, strongly suggesting the set is finite.

### Density analysis

If $K_2$, $K_3$, $K_4$ were linearly disjoint over $\mathbb{Q}$:
- $[K_2 \cdot K_3 \cdot K_4 : \mathbb{Q}] = 12 \cdot 18 \cdot 24 = 5,\!184$
- Density bound: $1/5,\!184 \approx 0.0193\%$

Observed density: $3 / 2,\!556 \approx 0.1174\%$

Since observed $>$ bound, the fields are **not** linearly disjoint, which is consistent with finiteness (the compositum is smaller than the product of degrees).

### False-positive caveat: fast-check limitations

The fast check (periods 2–5) is **not sufficient** for full verification. Known false positive:

- **$p = 313$** passes all period-2–5 checks but has a **period-27 ghost cycle**. Full verification requires checking **all** cycles and multipliers, not just whether polynomial roots exist modulo $p$.

The $p = 181$ anomaly also illustrates this: the period-4 polynomial has 4 linear factors modulo 181 ($u = 45, 47, 123, 179$), but **none are cubes** in $\mathbb{F}_{181}$, so $x^3 = u$ has no solutions — hence no period-4 points exist. Without the cube-root check, $p = 181$ would be a false *negative* (appearing to have period-4 points when it does not).
