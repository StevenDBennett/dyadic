"""
03_weight_stability.py — Demonstrate WeightStabilityDiagnostics

Usage: python examples/03_weight_stability.py
"""

import numpy as np
from dyadic_math import WeightStabilityDiagnostics


def main():
    np.random.seed(42)

    print("=" * 60)
    print("WeightStabilityDiagnostics — 2-Adic Weight Stability")
    print("=" * 60)

    W = np.random.randint(0, 256, size=(8, 8), dtype=np.int64)
    thermo = WeightStabilityDiagnostics(k=8)

    result = thermo.analyse(W)
    print("  Weight analysis:")
    for k, v in result.items():
        print(f"    {k}: {v:.4f}")

    print()

    print("=" * 60)
    print("Precision Sweep Analysis")
    print("=" * 60)

    weights = np.random.randint(1, 256, size=20, dtype=np.int64)
    thermo = WeightStabilityDiagnostics.from_precision_sweep(weights, k_range=range(4, 13), k=8)
    thermo.compute()

    print(f"  Stable weights: {thermo.stable_weights()}")
    print(f"  Ghost-affected weights: {thermo.ghost_weights()}")

    hist = thermo.cliff_histogram()
    if hist:
        print("  Cliff histogram:")
        for k, count in sorted(hist.items()):
            bar = "#" * count
            print(f"    k*={k:>2}: {count:>2} {bar}")

    print()
    print(thermo.report())

    print()

    print("=" * 60)
    print("Mersenne Cliff Score")
    print("=" * 60)

    for w in [3, 7, 15, 31, 63, 127, 255]:
        score = thermo.mersenne_cliff_score(w)
        print(f"  w={w:>4}: Mersenne cliff score = {score}")

    print()
    print("Done.")


if __name__ == "__main__":
    main()
