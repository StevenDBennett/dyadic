"""
crt.py
------
Chinese Remainder Theorem dual system:  Z/(2^k · p)Z  for odd prime p.

Combines the 2-adic dual-view decomposition with a modular (prime-field)
residue, allowing joint analysis of weights in the product ring.
"""

from __future__ import annotations

import numpy as np
from dyadic_core import DualNumber, modinv_newton, valuation


def _primitive_root(p: int) -> int | None:
    """Find the smallest primitive root modulo prime p."""
    if p == 2:
        return 1
    phi = p - 1
    factors = set()
    n = phi
    d = 2
    while d * d <= n:
        if n % d == 0:
            factors.add(d)
            while n % d == 0:
                n //= d
        d += 1
    if n > 1:
        factors.add(n)

    for g in range(2, p):
        ok = True
        for q in factors:
            if pow(g, phi // q, p) == 1:
                ok = False
                break
        if ok:
            return g
    return None


def _prime_dlog(a: int, p: int, g_p: int) -> int:
    """Brute-force discrete logarithm modulo prime p."""
    cur = 1
    for e in range(p):
        if cur == a:
            return e
        cur = (cur * g_p) % p
    return 0


class CRTDualNumber:
    """
    Element of Z/(2^k · p)Z with dual-coordinate representation.

    Attributes
    ----------
    component_2 : DualNumber — the 2-adic part
    residue_p  : int — residue modulo p
    dlog_p     : int — discrete-log of residue_p modulo p (base g_p)
    """

    def __init__(self, n: int, k: int, p: int, g_p: int) -> None:
        self.k = k
        self.p = p
        self.g_p = g_p
        self.mod_full = (1 << k) * p
        self.component_2 = DualNumber(n % (1 << k), k)
        self.residue_p = n % p
        self.dlog_p = _prime_dlog(self.residue_p, p, g_p) if self.residue_p != 0 else -1

    def verify(self) -> bool:
        """CRT round-trip check."""
        return self.component_2.verify()

    def __repr__(self) -> str:
        c2 = self.component_2
        if c2.is_zero:
            return f"CRTDualNumber(0, k={self.k}, p={self.p})"
        sign = "-" if c2.alpha else "+"
        return (
            f"CRTDualNumber(k={self.k}, p={self.p})"
            f"  =  2^{c2.v} · {sign}5^{c2.e}  (mod 2^{self.k})"
            f"  ×  g_p^{self.dlog_p}  (mod {self.p})"
        )


class CRTDualProcessor:
    """
    Arithmetic processor for CRTDualNumbers over Z/(2^k · p)Z.
    """

    def __init__(self, k: int, p: int, g_p: int | None = None) -> None:
        self.k = k
        self.p = p
        self.mod2 = 1 << k
        self.mod_full = self.mod2 * p
        if g_p is None:
            g_p = _primitive_root(p)
            if g_p is None:
                raise ValueError(f"No primitive root found for p={p}")
        self.g_p = g_p

    def crt_reconstruct(self, r2: int, rp: int) -> int:
        """
        CRT: combine residue r2 (mod 2^k) and rp (mod p) into
        a single integer mod 2^k·p.
        """
        inv_p = modinv_newton(self.p, self.k)
        t = ((r2 - rp) * inv_p) & (self.mod2 - 1)
        return (t * self.p + rp) % self.mod_full

    def mul(self, a_val: CRTDualNumber, b_val: CRTDualNumber) -> CRTDualNumber:
        """Multiply two CRTDualNumbers directly."""
        n2 = (a_val.component_2.value * b_val.component_2.value) % self.mod2
        np_val = (a_val.residue_p * b_val.residue_p) % self.p
        n_full = self.crt_reconstruct(n2, np_val)
        return CRTDualNumber(n_full, self.k, self.p, self.g_p)

    def product(self, weights: list[int]) -> CRTDualNumber:
        """Product of raw integers in CRT space."""
        nums = [CRTDualNumber(w, self.k, self.p, self.g_p) for w in weights]
        return self.cycle_product(nums)

    def cycle_product(self, numbers: list[CRTDualNumber]) -> CRTDualNumber:
        """Product of a list of CRTDualNumbers."""
        if not numbers:
            return CRTDualNumber(1, self.k, self.p, self.g_p)
        prod = numbers[0]
        for num in numbers[1:]:
            prod = self.mul(prod, num)
        return prod

    def convergence_ratio_2adic(self, product: CRTDualNumber) -> float:
        """Convergence ratio of the 2-adic component of product."""
        from .basin import BasinExplorer

        a_val = product.component_2.value
        if a_val == 0 or (a_val & 1) == 0:
            return 0.0
        try:
            explorer = BasinExplorer(self.k, 5, a_val)
            portrait = explorer.portrait()
            n_converged = len(portrait["converged"])
            total = n_converged + len(portrait["cycle"])
            return n_converged / total if total > 0 else 0.0
        except Exception:
            return 0.0

    def __repr__(self) -> str:
        return f"CRTDualProcessor(k={self.k}, p={self.p}, g_p={self.g_p})"


def combined_stability(
    k: int, p: int, num_cycles: int = 50, cycle_length: int = 4
) -> dict[str, float]:
    """
    Randomised test: correlation between 2-adic depth of a CRT cycle product
    and the 2-adic valuation of the change under an additive 2ᵗ perturbation.

    Unlike the original bit-flip approach (which always gave v₂(Δ)=0 because
    XOR by a power of 2 always produces an odd difference), this uses an
    additive perturbation w → w + 2ᵗ, producing a graded v₂(Δ) that varies
    meaningfully.  The stability measure is the v₂ of the original product's
    2-adic component — deeper in the filtration → more stable.

    Returns dict with pearson_r, n_samples, and mean_v2.
    """
    proc = CRTDualProcessor(k, p)
    v2_orig_vals: list[float] = []
    delta_v2s: list[float] = []

    for _ in range(num_cycles):
        weights = [np.random.randint(0, proc.mod_full) for _ in range(cycle_length)]
        product = proc.cycle_product([CRTDualNumber(w, k, p, proc.g_p) for w in weights])
        v2_orig_val = valuation(product.component_2.value)
        if v2_orig_val is None:
            v2_orig_val = 0

        # Additive 2-adic perturbation: add 2^t to one weight
        idx = np.random.randint(0, cycle_length)
        t = np.random.randint(0, k)
        perturbed = (weights[idx] + (1 << t)) % proc.mod_full
        weights_pert = weights[:]
        weights_pert[idx] = perturbed
        product_pert = proc.cycle_product([CRTDualNumber(w, k, p, proc.g_p) for w in weights_pert])

        delta = abs(product.component_2.value - product_pert.component_2.value)
        v2_delta_val: int = 0
        if delta > 0:
            v2_d = valuation(delta)
            v2_delta_val = v2_d if v2_d is not None else 0

        v2_orig_vals.append(float(v2_orig_val))
        delta_v2s.append(float(v2_delta_val))

    v2_arr = np.array(v2_orig_vals)
    deltas_arr = np.array(delta_v2s)
    std_v2 = v2_arr.std()
    std_d = deltas_arr.std()

    if std_v2 > 0 and std_d > 0:
        pearson_r = float(np.corrcoef(v2_arr, deltas_arr)[0, 1])
    else:
        pearson_r = 0.0

    return {
        "pearson_r": pearson_r,
        "n_samples": num_cycles,
        "mean_v2": float(v2_arr.mean()),
    }
