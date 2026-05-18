"""
clean_primes.py
---------------
Clean-prime detection for p-adic Newton dynamics.

A prime p is *clean* for N(x) = (2x^3+1)/(3x^2) if it admits no
periodic points of period 2, 3, or 4 with multiplier μ ≠ 1 (mod p).

The known clean primes are {7, 103, 181} — conjectured to be complete.
"""
from __future__ import annotations


__all__ = ["is_cube", "tonelli_shanks", "check_quadratic_cube_roots"]


def is_cube(a: int, p: int) -> bool:
    """Check whether a is a cube modulo the prime p (p ≡ 1 mod 3)."""
    if a % p == 0:
        return True
    return pow(a % p, (p - 1) // 3, p) == 1


def tonelli_shanks(n: int, p: int) -> int | None:
    """
    Tonelli–Shanks modular square root.

    Returns x such that x² ≡ n (mod p), or None if n is not a quadratic
    residue modulo the odd prime p.
    """
    if pow(n, (p - 1) // 2, p) != 1:
        return None
    if p % 4 == 3:
        return pow(n, (p + 1) // 4, p)
    # p ≡ 1 (mod 4)
    q, s = p - 1, 0
    while q % 2 == 0:
        q //= 2
        s += 1
    z = 2
    while pow(z, (p - 1) // 2, p) != p - 1:
        z += 1
    c = pow(z, q, p)
    x = pow(n, (q + 1) // 2, p)
    t = pow(n, q, p)
    m = s
    while t != 1:
        i = 1
        while i < m and pow(t, 2 ** i, p) != 1:
            i += 1
        b = pow(c, 2 ** (m - i - 1), p)
        x = (x * b) % p
        t = (t * b * b) % p
        c = (b * b) % p
        m = i
    return x


def check_quadratic_cube_roots(a: int, b: int, c: int, p: int) -> bool:
    """
    Check whether the quadratic ax² + bx + c ≡ 0 (mod p) has a root
    that is a cube modulo p (i.e. lies in F_p³).

    Returns True if at least one root is a cube modulo p.
    """
    a, b, c = a % p, b % p, c % p
    if a == 0:
        if b == 0:
            return False
        return is_cube((-c * pow(b, -1, p)) % p, p)
    disc = (b * b - 4 * a * c) % p
    if pow(disc, (p - 1) // 2, p) != 1:
        return False
    sqrt_d = tonelli_shanks(disc, p)
    if sqrt_d is None:
        return False
    inv = pow(2 * a, -1, p)
    r1 = ((-b + sqrt_d) * inv) % p
    r2 = ((-b - sqrt_d) * inv) % p
    return is_cube(r1, p) or is_cube(r2, p)
