"""
01_dual_number_basics.py — Demonstrate DualNumber, TwoAdicProcessor, modinv

Usage: python examples/01_dual_number_basics.py
"""

from dyadic_core import DualNumber, TwoAdicProcessor, modinv_newton, padic_exp, padic_log


def main():
    print("=" * 60)
    print("DualNumber — 2-adic Decomposition")
    print("=" * 60)

    for n in [0, 1, 3, 5, 7, 42, 255, 256]:
        dn = DualNumber(n, k=16)
        status = "✓" if dn.verify() else "✗"
        print(f"  {status} DualNumber({n:>4}, k=16): v={dn.v:>2}, α={dn.alpha}, e={dn.e:>5}")

    print()

    print("=" * 60)
    print("TwoAdicProcessor — Coordinate-Space Arithmetic")
    print("=" * 60)

    proc = TwoAdicProcessor(k=16)
    a = DualNumber(3, k=16)
    b = DualNumber(5, k=16)
    c = proc.mul(a, b)
    print(f"  3 × 5 = {c.value}  (expected 15)")

    inv3 = proc.inv(a)
    check = proc.mul(a, inv3).value
    print(f"  3 × 3⁻¹ = {check}  (expected 1)")

    neg3 = proc.pow(a, -1)
    print(f"  3⁻¹ (via pow) = {neg3.value}")

    print()

    print("=" * 60)
    print("modinv_newton — Modular Inverse")
    print("=" * 60)

    for a in [3, 5, 7, 9, 11]:
        inv = modinv_newton(a, k=16)
        mask = (1 << 16) - 1
        check = (a * inv) & mask
        status = "✓" if check == 1 else "✗"
        print(f"  {status} {a}⁻¹ mod 2¹⁶ = {inv:>6}  (check: {a}×{inv} mod = {check})")

    print()

    print("=" * 60)
    print("2-Adic Exponential and Logarithm")
    print("=" * 60)

    for x in [4, 8, 12]:
        exp_x = padic_exp(x, k=16)
        log_exp = padic_log(exp_x, k=16)
        print(f"  exp({x:>2}) = {exp_x:>6},  log(exp) = {log_exp}")

    print()
    print("All checks passed.")


if __name__ == "__main__":
    main()
