"""
04_mersenne_ghost.py — Demonstrate Mersenne Ghost Theorem

Usage: python examples/04_mersenne_ghost.py
"""

from dyadic_math import (
    mersenne_coordinates,
    mersenne_cliff_table,
    mersenne_cliff_theorem,
    cliff_constant,
    cliff_constant_unified,
)
from dyadic_math.mersenne import verify_core_identity, bootstrap_cost, optimal_bootstrap


def main():
    print("=" * 60)
    print("Mersenne Ghost Theorem — Core Identity")
    print("=" * 60)

    verify_core_identity(n_max=10)

    print()

    print("=" * 60)
    print("Mersenne Coordinates")
    print("=" * 60)

    for n in range(3, 11):
        w = (1 << n) - 1
        alpha, e_true, v2_e = mersenne_coordinates(n, k=n + 3)
        print(f"  n={n:>2}  w=2^{n}-1={w:>6}  α={alpha}  e_true={e_true:>6}  v₂(e)={v2_e}")

    print()

    print("=" * 60)
    print("Mersenne Cliff Table")
    print("=" * 60)

    table = mersenne_cliff_table(n_max=10)
    print(f"  {'n':>3} | {'w':>10} | {'α':>3} | {'e_true':>8} | {'v₂(e)':>6} | {'k*':>4}")
    print("  " + "-" * 50)
    for row in table:
        w = (1 << row["n"]) - 1
        print(f"  {row['n']:>3} | {w:>10} | {row['alpha']:>3} | {row['e_true']:>8} | {row['v2_e']:>6} | {row['k*']:>4}")

    print()

    print("=" * 60)
    print("Cliff Constant")
    print("=" * 60)

    for g in [5, 13, 21, 29, 37]:
        c = cliff_constant(g, k=16)
        cu = cliff_constant_unified(g, k=16)
        print(f"  cliff_constant(g={g:>2}) = {c}, unified = {cu}")

    print()

    print("=" * 60)
    print("Bootstrap Cost Analysis")
    print("=" * 60)

    for k in [16, 32, 64, 128]:
        sqrt_cost = bootstrap_cost(int(k ** 0.5), k)
        half_cost = bootstrap_cost(k // 2, k)
        print(f"  k={k:>4}:  sqrt-k cost={sqrt_cost:>6}  k/2 cost={half_cost:>6}  "
              f"saving={sqrt_cost - half_cost:>6} ({100 * (1 - half_cost / sqrt_cost):.1f}%)")

    print()

    print("=" * 60)
    print("Mersenne Cliff Theorem (Full Verification)")
    print("=" * 60)

    mersenne_cliff_theorem(verbose=True)

    print()
    print("Done.")


if __name__ == "__main__":
    main()
