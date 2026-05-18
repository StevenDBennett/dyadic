"""
02_basin_explorer.py — Demonstrate BasinExplorer ghost detection

Usage: python examples/02_basin_explorer.py
"""

from dyadic_math import BasinExplorer, mersenne_cliff_table


def main():
    print("=" * 60)
    print("Basin Explorer — Newton Basin Landscape")
    print("=" * 60)

    for target in [1, 3, 5, 7]:
        explorer = BasinExplorer(k=8, g=5, target_a=target)
        portrait = explorer.portrait()
        n_conv = len(portrait["converged"])
        n_cycle = len(portrait["cycle"])
        n_div = len(portrait["diverged"])
        fate = "GHOST" if n_cycle > 0 else "CLEAN"
        print(f"  target={target:>2} | α={explorer.alpha} | "
              f"converged={n_conv:>3} cycle={n_cycle:>3} diverged={n_div:>3} | {fate}")

    print()

    print("=" * 60)
    print("Basin Classification — Single Seed")
    print("=" * 60)

    for target, e0 in [(1, 0), (3, 0), (3, 16), (7, 32)]:
        explorer = BasinExplorer(k=8, g=5, target_a=target)
        fate, value, path = explorer.classify(e0)
        print(f"  target={target:>2}, seed e₀={e0:>2}: fate={fate:>12}, value={value}")

    print()

    print("=" * 60)
    print("Mersenne Cliff Table")
    print("=" * 60)

    table = mersenne_cliff_table(n_max=8)
    print(f"  {'n':>3} | {'w':>10} | {'α':>3} | {'e_true':>8} | {'v₂(e)':>6} | {'k*':>4}")
    print("  " + "-" * 50)
    for row in table:
        w = (1 << row["n"]) - 1
        print(f"  {row['n']:>3} | {w:>10} | {row['alpha']:>3} | {row['e_true']:>8} | {row['v2_e']:>6} | {row['k*']:>4}")

    print()
    print("Done.")


if __name__ == "__main__":
    main()
