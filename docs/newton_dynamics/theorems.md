> This document is part of the dyadic (3-package monorepo) project, ported from dual-view.

# Theorems

Four theorems have been proven for the Newton map $N(x) = (2x^3+1)/(3x^2)$.

---

## Theorem 1: $\Phi_n^*(0) = 2^{\mu_3(n)/6}$ (Proven)

**Proof:**

The Newton iterates $N^d(x) = A_d(x)/B_d(x)$ satisfy:
$$A_{d+1} = 2A_d^3 + B_d^3, \qquad B_{d+1} = 3A_d^2 B_d$$
with $A_0 = x$, $B_0 = 1$.

At $x = 0$: $A_d(0) = 2^{(3^{d-1}-1)/2}$ for $d \ge 1$, and $B_d(0) = 0$.

The dynatomic polynomial at $x = 0$ is:
$$\Phi_n^*(0) = \prod_{d|n} A_d(0)^{\mu(n/d)} = 2^{\sum_{d|n} \frac{3^{d-1}-1}{2} \cdot \mu(n/d)}$$

The exponent simplifies:
$$\sum_{d|n} \frac{3^{d-1}-1}{2} \cdot \mu(n/d) = \frac{1}{2}\sum 3^{d-1}\mu(n/d) - \frac{1}{2}\sum\mu(n/d) = \frac{1}{6}\mu_3(n) - 0$$

since $\sum_{d|n} \mu(n/d) = 0$ for $n > 1$. Therefore $\Phi_n^*(0) = 2^{\mu_3(n)/6}$. ∎

---

## Theorem 2: $\Phi_n^*(1) = 3^{\mu_3(n)/2}$ (Proven)

**Proof (Lemma 1 — derivative identity):**

Define $F_d(x) = A_d(x) - x \cdot B_d(x)$. Then $N^d(x) - x = F_d(x)/B_d(x)$.

$F_d(1) = 0$ since $N^d(1) = 1$. The derivatives satisfy:
$$F_d'(x) = A_d'(x) - B_d(x) - x \cdot B_d'(x)$$

From the recurrence:
$$A_{d+1}' = 6A_d^2 A_d' + 3B_d^2 B_d'$$
$$B_{d+1}' = 6A_d A_d' B_d + 3A_d^2 B_d'$$

At $x=1$: $A_d(1) = B_d(1)$, so $A_{d+1}'(1) = B_{d+1}'(1)$. By induction with base case $A_1'(1) = B_1'(1) = 6$:
$$A_d'(1) = B_d'(1) \text{ for all } d$$

Therefore $F_d'(1) = B_d'(1) - B_d(1) - B_d'(1) = -B_d(1)$. ∎

**Proof (Lemma 2 — B_d(1) closed form):**

$B_d(1) = 3^{(3^d-1)/2}$ for all $d \ge 1$.

*Proof by induction*:
- Base $d=1$: $B_1(x) = 3x^2$, so $B_1(1) = 3 = 3^{(3-1)/2}$.
- Step: $B_{d+1}(1) = 3 \cdot A_d(1)^2 \cdot B_d(1) = 3 \cdot B_d(1)^3 = 3 \cdot (3^{(3^d-1)/2})^3 = 3^{1 + (3^{d+1}-3)/2} = 3^{(3^{d+1}-1)/2}$. ∎

**Final computation:**

Near $x=1$, each $F_d$ has a simple zero. By L'Hôpital:
$$\lim_{x\to 1} \frac{F_d(x)}{F_e(x)} = \frac{F_d'(1)}{F_e'(1)} = \frac{B_d(1)}{B_e(1)}$$

After clearing the denominator factors (which cancel for $e < n$):
$$\Phi_n^*(1) = \prod_{d|n} B_d(1)^{\mu(n/d)} = \prod_{d|n} \left(3^{(3^d-1)/2}\right)^{\mu(n/d)}$$
$$= 3^{\frac{1}{2}\sum 3^d\mu(n/d) - \frac{1}{2}\sum\mu(n/d)} = 3^{\mu_3(n)/2}$$

since $\sum_{d|n} \mu(n/d) = 0$ for $n > 1$. ∎

---

## Theorem 3: $\prod\mu = 6^{\mu_3(n)/6}$ (Universal Product Formula — Proven)

**Proof:**

For a period-$n$ cycle $\{x_i\}$ with $u_i = x_i^3$:
$$\mu = \prod_{i=0}^{n-1} N'(x_i) = \left(\frac{2}{3}\right)^n \prod_{i=0}^{n-1} \frac{u_i - 1}{u_i}$$

For all cycles combined ($c = \mu_3(n)/(3n)$ cycles):
$$\prod\mu = \left(\frac{2}{3}\right)^{n \cdot c} \cdot \prod_{\text{all roots}} \frac{u_i - 1}{u_i} = \left(\frac{2}{3}\right)^{\mu_3(n)/3} \cdot \frac{\Phi_n^*(1)}{\Phi_n^*(0)}$$

Using Theorems 1 and 2:
$$\frac{\Phi_n^*(1)}{\Phi_n^*(0)} = \frac{3^{\mu_3(n)/2}}{2^{\mu_3(n)/6}} = \left(\frac{27}{2}\right)^{\mu_3(n)/6}$$

Therefore:
$$\prod\mu = \left(\frac{2}{3}\right)^{\mu_3(n)/3} \cdot \left(\frac{27}{2}\right)^{\mu_3(n)/6} = 6^{\mu_3(n)/6}$$

| Period | $\mu_3(n)$ | $\prod\mu$ |
|--------|------------|------------|
| 2 | 6 | $6^1 = 6$ |
| 4 | 72 | $6^{12} = 2,\!176,\!782,\!336$ |
| 5 | 240 | $6^{40}$ |

∎

---

## Theorem 4: $\mu_u = \mu_x$ (Period-4 Multiplier Identity — Proven)

**Proof:**

For a 4-cycle $\{u_i\}$ under $M(u) = \frac{(2u+1)^3}{27u^2}$, the cycle equation gives:
$$\prod_{i=0}^3 (2u_i + 1) = 81 \prod_{i=0}^3 u_i$$

The multiplier in $u$-dynamics is $\mu_u = \frac{16}{81} \prod (u_i - 1)/u_i$. The multiplier in $x$-dynamics is $\mu_x = \left(\frac{2}{3}\right)^4 \prod (u_i - 1)/u_i = \frac{16}{81} \prod (u_i - 1)/u_i$. Therefore $\mu_u = \mu_x$. ∎

---

## Bonus: 3-Adic Quadratic Convergence

The Newton map $N(x) = (2x^3+1)/(3x^2)$ has superattracting fixed points at the cube roots of unity ($N'(\zeta) = 0$ for $\zeta^3 = 1$). Near $x = 1$:

$$N(1+\varepsilon) = \frac{2(1+\varepsilon)^3+1}{3(1+\varepsilon)^2} \approx 1 + \varepsilon^2$$

so $v_3(N(1+\varepsilon)-1) = 2 \cdot v_3(\varepsilon)$ — the map is **quadratically convergent** in the 3-adic metric. This is consistent with LTE:

$$v_3(4^n - 1) = 1 + v_3(n) \quad (n \text{ divisible by } 3)$$

and the general theory of superattracting fixed points.
